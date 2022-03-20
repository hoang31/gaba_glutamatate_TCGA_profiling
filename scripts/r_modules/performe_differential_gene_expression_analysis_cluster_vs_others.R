
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

## load ensembl id
ensembl_id <- fread(
    file = snakemake@params[["ensembl_id"]],
    sep = ","
)

## load external scripts
source(snakemake@params[["utils_R"]])

## create the output directory
#dir.create(snakemake@output[["DEgene_results_dir"]])


###################################################
###################################################


##### REFORMATE THE DATA


## rename the column containing the gene name
colnames(d_expression)[colnames(d_expression) == "genes"] = "Gene_Name"

## remove useless rows
row_to_remove <- c("__alignment_not_unique", "__ambiguous", "__no_feature")
d_expression <- d_expression[!(Gene_Name %in% row_to_remove), ]

## create the directory that will contain the data
dir.create(snakemake@output[["results_TCGA_IDHall_cluster_vs_others_dir"]]) 

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
    CLUSTER_GROUP_INPUT, # data table containing two columns : sample_id and cluster
    cluster_of_interest # name of the cluster that have to be compared with the others
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

    ## take only the samples of the cluster "MOXTE"
    coldata_output[cluster %in% c(cluster_of_interest), cluster := cluster_of_interest]
    coldata_output[!(cluster %in% c(cluster_of_interest)), cluster := "OTHERS"]

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
    expression_matrix_input[
        ,
        lapply(.SD, function(x) round(x)),
        .SDcols = colnames(expression_matrix_input)[colnames(expression_matrix_input) != "Gene_Name"]
    ]

    setDF(
        expression_matrix_input,
        rownames = unlist(expression_matrix_input[, "Gene_Name"]
    ))

    ## remove the Gene_Name column 
    expression_matrix_input[, "Gene_Name"] = NULL

    ## return the expression matrix
    return(expression_matrix_input)
}


###################################################
###################################################

##### DIFFERENTIAL GENE EXPRESSION ANALYSIS

## function for apply DESEQ2 normalization
deseq2_analysis <- function(
    CTS_INPUT, # corresponding to the expression data
    COLDATA_INPUT,
    cluster_name, # name of the cluster of interest
    file_path # corresponding to the path where the file will be written
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

    ## differential expression analysis
    dds <- DESeq(
        dds,
        parallel = FALSE
    )

    ## Extract the genes that are differentially expressend between mixed and others clusters
    res <- results(
        dds,
        contrast = c('cluster', cluster_name, 'OTHERS'),
        independentFiltering = F,
        pAdjustMethod = "bonferroni",
    )

    ## transform to data frame
    res <- as.data.table(as.data.frame(res), keep.rownames = 'Gene_Name')[order(padj, decreasing = F),]

    ## transform the ensembl id to gene symbols
    res[, Gene_Name := ensemblID_to_geneSymbole(Gene_Name, ensembl_id)]

    ## write the output results
    fwrite(
        res,
        file_path,
        sep = ","
    )
    
}


###################################################
###################################################

##### CALL THE FUNCTION

## retreive the cluster names into a vector
cluster_name_vector <- unique(unlist(d_clusters[, cluster]))



## do a loop for all the cluster names
for (i_cluster in seq(1, length(cluster_name_vector), 1)) {
    
    ## create the col data
    coldata <- create_coldata(
        EXPRESSION_MATRIX_INPUT = d_expression,
        CLUSTER_GROUP_INPUT = d_clusters,
        cluster_of_interest = cluster_name_vector[i_cluster]
    )

    ## retreinve the expression data
    d_expression_filtered <- extract_expression_data(
        EXPRESSION_MATRIX_INPUT = d_expression,
        COL_DATA_INPUT = coldata
    )

    ## generate the path of the file to be written
    file_path <- paste(
        snakemake@output[["results_TCGA_IDHall_cluster_vs_others_dir"]],
        "/",
        "results_",
        cluster_name_vector[i_cluster],
        "_vs_others",
        ".csv",
        sep = ""
    )

    ## perform the differential gene expression analysis 
    deseq2_output <- deseq2_analysis(
        CTS_INPUT = d_expression_filtered,
        COLDATA_INPUT = coldata,
        cluster_name = cluster_name_vector[i_cluster],
        file_path = file_path
    )

}

#break
### create the col data
#coldata <- create_coldata(
#    EXPRESSION_MATRIX_INPUT = d_expression,
#    CLUSTER_GROUP_INPUT = d_clusters,
#    cluster_of_interest = "NT-1"
#)

### retreinve the expression data
#d_expression_filtered <- extract_expression_data(
#    EXPRESSION_MATRIX_INPUT = d_expression,
#    COL_DATA_INPUT = coldata
#)

### perform the differential gene expression analysis 
#deseq2_output <- deseq2_analysis(
#    CTS_INPUT = d_expression_filtered,
#    COLDATA_INPUT = coldata
#)



###################################################
###################################################



