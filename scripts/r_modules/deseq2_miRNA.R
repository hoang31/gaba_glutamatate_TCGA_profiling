
###########################################
## Do the differential miRNA expression using DEseq2
###########################################


###########################################
## Load the libraries
###########################################


library(data.table)
library(tidyverse)
library(BiocParallel)
suppressMessages(library(DESeq2))
register(MulticoreParam(detectCores() - 2))


###########################################
## Load the data
###########################################


## load the rawcounts of the miRNA expression data
miRNA_expression_data_dt <- fread(
    snakemake@input[["miRNA_expression_data_rawcount"]],
    sep = ","
)

## load the cluster data
cluster_data_dt <- fread(
    snakemake@input[["cluster_group"]],
    sep = ","

)

## load the utils script
source(snakemake@params[["utils"]])

## create the directory that will contain the differentially miRNA expression analysis results
dir.create(snakemake@output[["results_deseq2_miRNA"]])


#####
## read the healthy sample mir expression data
#####

## read the path
healty_path <- list.files(
    snakemake@input[["mir_expression_healthy_dir"]],
    full.names = TRUE
)

## read the table for each path
healthy_mir_expression_dt <- lapply(
    healty_path,
    function(x) {

        ## get the name of the sample
        sample_name_value <- str_match(
            x,
            pattern = "data/expression_data/TCGA_mir_healthy_samples/(.*?).txt"
        )[1,2]

        ## read the table
        output <- fread(
            x,
            sep = "\t"
        )[, c('miRNA_ID', 'read_count'), with = F][, setnames(.SD, c("miRNA_ID", "read_count"), c("gene_name", sample_name_value))]

        return(output)
    }
)

## merging the data together
healthy_mir_expression_dt <- Reduce(
    x = healthy_mir_expression_dt,
    f = merge
)

## adding the healthy sample id to the cluster data
cluster_data_dt <- rbind(
    cluster_data_dt,
    data.table(
        sample_id = colnames(healthy_mir_expression_dt)[colnames(healthy_mir_expression_dt) != "gene_name"],
        cluster = "Healthy"
    )
)


#####
## Add the healthy sample into the whole expression data
#####

print(dim(miRNA_expression_data_dt))

miRNA_expression_data_dt <- merge(
    miRNA_expression_data_dt, 
    healthy_mir_expression_dt,
    by = "gene_name",
    sort = F
)



###########################################
## Preparation of the data
###########################################


## do the transposition of the miRNA_expression_data_dt
miRNA_expression_data_dt <- transpose_datatable(
    miRNA_expression_data_dt,
    column_name = "gene_name",
    new_name = "sample_id"
)

## merge with the cluster data
miRNA_expression_data_dt <- merge(
    cluster_data_dt,
    miRNA_expression_data_dt,
    by = 'sample_id'
)


##########
## Prepare the information sample data
##########


## create the information data table from the miRNA expression data table
information_data_dt <- miRNA_expression_data_dt[,c('sample_id', 'cluster'), with = F]


## transform the data table to data frame
setDF(
    information_data_dt,
    rownames = unlist(miRNA_expression_data_dt[, 'sample_id', with = F])
)


##########
## Prepare the miRNA expression data
##########


## function for preparing the expression data from the information data
extract_expression_data <- function(
    expression_data_input,
    information_data_input
) {
    ## copy the expression matrix input
    expression_matrix_input <- copy(expression_data_input[,!c("cluster")])

    ### extract the exprssion related to the samples in the coldata
    #expression_matrix_input <- expression_matrix_input[, append("sample_id", rownames(information_data_input)), with = F]

    ## do the transposition
    expression_matrix_input <- transpose_datatable(
        expression_matrix_input,
        "sample_id",
        "gene_name"
    )

    ## transform all the value to numeric values
    expression_matrix_input_numeric <- expression_matrix_input[
        ,
        lapply(.SD, function(x) (as.numeric(x))),
        .SDcols = colnames(expression_matrix_input)[colnames(expression_matrix_input) != "gene_name"]
    ]

    setDF(
        expression_matrix_input_numeric,
        rownames = unlist(expression_matrix_input[, "gene_name"])
    )

    ## remove the gene_name column 
    expression_matrix_input_numeric[, "gene_name"] = NULL

    ## return the expression matrix
    return(expression_matrix_input_numeric)
}

## call the function for preparing the data
miRNA_expression_data_formated <- extract_expression_data(
    miRNA_expression_data_dt,
    information_data_dt
)


###########################################
## Call DESeq2
###########################################


##
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

    ## differential expression analysis (parallelization)
    dds <- DESeq(
        dds,
        parallel = TRUE
    )

    ## estimate the size factor for the normalization
    size_factor <- estimateSizeFactors(dds)
    ## variance stabilizing transformation
    vsd <- assay(varianceStabilizingTransformation(size_factor, blind = FALSE, fitType = "parametric"))

    ## extract the normalized counts
    normalized_cts <- data.table(vsd, keep.rownames = "gene_name")

    ## write the normalized data into a file
    fwrite(
        normalized_cts,
        snakemake@output[["deseq2_miRNA_normalized_expression"]],
        sep = ","
    )

    ## take the label name 
    cluster_labels <- unique(unlist(COLDATA_INPUT[, "cluster"]))




    ## remomving the healthy condition from the DE gene expression analysis
    cluster_labels <- cluster_labels[cluster_labels != "Healthy"]

    print(cluster_labels)






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

        ## create the file output name
        file_output_name <- paste(
            paste(
                combination,
                sep = "_"
            ),
            "csv",
            sep = "."
        )
 
        ## create the file output path
        file_output_path <- paste(  
            snakemake@output[["results_deseq2_miRNA"]],
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
            "gene_name"
        )][
            order(padj, decreasing = F),
        ]

        ## write the output results
        fwrite(
            res,
            file_output_path,
            sep = ","
        )

    }
}

## call the function for the deseq2 analysis
deseq2_analysis(
    miRNA_expression_data_formated,
    information_data_dt
)

