
###################################################
###################################################


## Differential gene expression analysis


###################################################
###################################################


## LOAD LIBRARIES
library(data.table)
library(tidyverse)
suppressMessages(library(DESeq2))


###################################################
###################################################


##### LOAD DATA

## load expression data
d_expression <- fread(
    file = snakemake@input[["TCGA_expression"]],
    sep = ",",
    header = T
)

## load clinical data
d_clinical<- fread(
    file = snakemake@input[["TCGA_clinical"]],
    sep = ",",
    header = T
)

## load the sample clusters
d_clusters <- fread(
    file = snakemake@input[["cluster_group"]],
    sep = ","
)

table(d_clusters[, "cluster"])

## retreive the clustering name
cluster_name <- str_match(
    snakemake@input[["cluster_group"]],
    pattern = "others/(.*?)_clusters.csv"
)[,2]


## load the genes of interest
gene_of_interest <- fread(
    file = snakemake@input[["gene_of_interest"]],
    sep = ","
)

## load ensembl id
ensembl_id <- fread(
    file = snakemake@params[["ensembl_id"]],
    sep = ","
)

## load external scripts
source(snakemake@params[["utils_R"]])

## create the output directory
dir.create(snakemake@output[["DEgene_results_dir"]])

###################################################
###################################################


##### REFORMATE THE DATA


## rename the column containing the gene name
colnames(d_expression)[colnames(d_expression) == "genes"] = "Gene_Name"

## remove useless rows
row_to_remove <- c("__alignment_not_unique", "__ambiguous", "__no_feature")
d_expression <- d_expression[!(Gene_Name %in% row_to_remove), ]


###################################################
###################################################


##### CREATE COLDATA


##### Function for formate the data for the deseq2 object

## function for create the "coldata" data frame containing all the information of the samples
## example of EXPRESSION_MATRIX_INPUT :
## Gene_Name    sample1   sample2   sample3
## gene1        x         x         x
## gene2        x         x         x
## gene2        x         x         x


create_coldata <- function(
    EXPRESSION_MATRIX_INPUT,
    CLUSTER_GROUP_INPUT # data table containing two columns : sample_id and cluster
) {

    ## initialize the coldata matrix
    coldata_output <- as.data.table(
        colnames(EXPRESSION_MATRIX_INPUT)[colnames(EXPRESSION_MATRIX_INPUT) != "Gene_Name"]
    )

    ## merge the two data tables
    coldata_output <- merge(
        coldata_output,
        CLUSTER_GROUP_INPUT,
        by.x = "V1",
        by.y = "sample_id",
        all.y = T,
        sort = F
    )



    ###################################################
    ###################################################
    ###################################################
    ################################################### 
    ###################################################
    ################################################### 





    # ## take only the samples of the cluster "MOXTE"
    # coldata_output[cluster %in% c("MIXTE", "IDH_MUT_codel", "IDH_MUT_noncodel"), cluster := "good_prognosis"]

    # coldata_output[cluster %in% c("IDH_WT_MIXTE", "IDH_WT"), cluster := "bad_prognosis"]





    ###################################################
    ################################################### 
    ###################################################
    ###################################################

    ## transform to data frame and rename the rownames with the sample id
    setDF(
        coldata_output,
        rownames = unlist(coldata_output[, "V1"]
    ))

    ## rename the column containing the sample_id
    colnames(coldata_output)[colnames(coldata_output) == "V1"] <- "sample_id" # rename the column

    ## transform to factor
    coldata_output$cluster <- factor(as.character(coldata_output$cluster))

    ## return the output
    return(coldata_output)
}


###################################################
###################################################


##### PREPARING THE EXPRESSION DATA


## function for extracting the expression data related to the samples in the coldata
extract_expression_data <- function(
    EXPRESSION_MATRIX_INPUT,
    COL_DATA_INPUT
) {
    ## copy the expression matrix input
    expression_matrix_input <- copy(EXPRESSION_MATRIX_INPUT)

    ## extract the exprssion related to the samples in the coldata
    expression_matrix_input <- expression_matrix_input[, append("Gene_Name", rownames(COL_DATA_INPUT)), with = F]

    ## transform all the value to numeric values
    expression_matrix_input_numeric <- expression_matrix_input[
        ,
        lapply(.SD, function(x) round(x)),
        .SDcols = colnames(expression_matrix_input)[colnames(expression_matrix_input) != "Gene_Name"]
    ]

    setDF(
        expression_matrix_input_numeric,
        rownames = unlist(expression_matrix_input[, "Gene_Name"]
    ))

    ## remove the Gene_Name column 
    expression_matrix_input_numeric[, "Gene_Name"] = NULL

    ## return the expression matrix
    return(expression_matrix_input_numeric)
}


###################################################
###################################################


##### DIFFERENTIAL GENE EXPRESSION ANALYSIS

## function for apply DESEQ2 normalization
deseq2_analysis <- function(
    CTS_INPUT, # corresponding to the expression data
    COLDATA_INPUT
) {
    ## create the dds object with the cts and coldata data
    dds <- DESeqDataSetFromMatrix(
        countData = CTS_INPUT,
        colData = COLDATA_INPUT,
        design = ~ cluster
    )

    ## keep only rows that have at least 10 reads total
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]



    # differential expression analysis (parallelization)
    # dds <- DESeq(
    #     dds,
    #     parallel = TRUE,
    #     BPPARAM=MulticoreParam(5)
    # )

    # differential expression analysis
    dds <- DESeq(
        dds,
        parallel = FALSE
    )

    ## take the label name 
    cluster_labels <- unique(unlist(COLDATA_INPUT[, "cluster"]))

    ## do all the 2 combination of the labels
    cluster_labels_combinations <- combn(
        x = cluster_labels,
        m = 2
    )

    ## for each combination
    for (combination_index in seq(1, ncol(cluster_labels_combinations), 1)) {

        cat("combination :", toString(cluster_labels_combinations[1, combination_index]), "vs", toString(cluster_labels_combinations[2, combination_index]), "\n")

        ## retreive the file path of the output based on the used combination

        ## create combination string
        combination <- paste(
            toString(cluster_labels_combinations[1, combination_index]),
            "vs",
            toString(cluster_labels_combinations[2, combination_index]),
            sep = "_"
        )
        print(combination)
        # break
        ## create the file output name
        file_output_name <- paste(
            paste(
                cluster_name,
                combination,
                sep = "_"
            ),
            "csv",
            sep = "."
        )
 
        ## create the file output path
        file_output_path <- paste(  
            snakemake@output[["DEgene_results_dir"]],
            file_output_name,
            sep = "/"
        )

        ## extract the results
        res <- results(
            dds,
            contrast = c("cluster", toString(cluster_labels_combinations[1, combination_index]), toString(cluster_labels_combinations[2, combination_index])),
            independentFiltering = F,
            pAdjustMethod = "bonferroni",
        )
        
        ## transform to data frame
        res <- as.data.frame(res)
        
        ## transforme to data table, rename the column containing the gene names and order the row by the padj
        res <- as.data.table(
            res,
            keep.rownames = T
        )[, setnames(
            .SD,
            'rn',
            "Gene_Name"
        )][
            order(padj, decreasing = F),
        ]

        ## write the output results
        fwrite(
            res,
            file_output_path,
            sep = ","
        )

         ## copy the results into a other data table
        res_filtered <- copy(res)   

        ## transform the ensembl id to gene symbols
        res_filtered[, Gene_Name := ensemblID_to_geneSymbole(Gene_Name, ensembl_id)]

        ## extract the res_filtered associated with the genes of interest
        res_filtered <- res_filtered[Gene_Name %in% unlist(gene_of_interest[,"gene_name"]), ]



        file_output_filtered_path <- paste(  
            snakemake@output[["DEgene_results_dir"]],
            '/filtered_',
            file_output_name,
            sep = ""
        )

        print(file_output_path)
        print(file_output_filtered_path)

        fwrite(
            res_filtered,
            file_output_filtered_path,
            sep = ","
        )
    }
}


###################################################
###################################################

##### CALL THE FUNCTION

## create the col data
coldata <- create_coldata(
    EXPRESSION_MATRIX_INPUT = d_expression,
    CLUSTER_GROUP_INPUT = d_clusters
)

## retreinve the expression data
d_expression_filtered <- extract_expression_data(
    EXPRESSION_MATRIX_INPUT = d_expression,
    COL_DATA_INPUT = coldata
)

## perform the differential gene expression analysis 
deseq2_output <- deseq2_analysis(
    CTS_INPUT = d_expression_filtered,
    COLDATA_INPUT = coldata
)

###################################################
###################################################



