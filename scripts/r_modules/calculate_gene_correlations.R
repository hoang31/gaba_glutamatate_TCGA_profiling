
##############################################
## Calculate the correlation between the genes of interest and the others genes
##############################################

##############################################
## Load the libraries
##############################################

library(data.table)
library(tidyverse)
library(foreach)


cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)


##############################################
## Load the data
##############################################

## load the genes of interest
gene_of_interest_vector <- scan(
    file = snakemake@input[["gene_of_interest"]],
    sep = "\n",
    what = character()
)

## load the table that contains the emsebl id and the gene symbol
ensembl_id <- fread(
    snakemake@params[['ensembl_id']]
)

## extract the ensembl id of the gene of interest
ensembl_id_gene_of_interest <- ensembl_id[gene_name %in% gene_of_interest_vector]

## load the expression data
expression_data <- fread(
    snakemake@input[["TCGA_expression_normalized"]]
)[, setnames(.SD, 'rn', 'gene_id')]

## load the directory that contains all the results of the differentially gene expression analysis
TCGA_deseq_res_dir <- snakemake@input[['TCGA_deseq_res_dir']]

##############################################
## ## extract the expression data of the general genes that are differentially expressed between the different group
##############################################

## list all the files contained by the deseq directory
DEgene_file_name <- list.files(TCGA_deseq_res_dir, full.names = T)

## remove the filtered DE genes
DEgene_file_name <- DEgene_file_name[-(grep(DEgene_file_name, pattern = "filtered"))]

## read all the results from the directory and put them into a same list
DEgene_data_list <- lapply(
    DEgene_file_name,
    function(x) {
        fread(
            x,
            sep = ","
        )
    }
)

## rename the data with the names
names(DEgene_data_list) <- DEgene_file_name

## rename the names of each sublist of data
names(DEgene_data_list) <- str_match(
    names(DEgene_data_list),
    pattern = paste(
        'IDHall',
        "_(.*?).csv",
        sep = ""
    )
)[,2]

#####
## Extract the DEgenes for each DESEQ results
#####

## initialize the p value and the logfoldchange thresolds
pval_threshold <- 0.05
logFC_threshold <- 2

## initialize the list that will contain the significant list
significant_gene_list <- list()

## for each subset of the merged_data, extract the genes associated with a pvalue lower than pval_threshold and put them into the significant_gene_list
for (i in seq(1, length(DEgene_data_list), 1)) {

    # print("###############################")

    ## extract the subset_data of the DEgene_data_list list
    subset_data <- as.data.table(DEgene_data_list[[i]])

    ## we have a problem when we are selecting the genes depending of the log fold change so I foound a workaround
    subset_data[, log2FoldChange := as.numeric(log2FoldChange)]
    subset_data <- (subset_data[(padj < pval_threshold),])
    subset_data <- subset_data[((log2FoldChange) > logFC_threshold), to_select := "yes"][order(log2FoldChange),]
    subset_data <- subset_data[((log2FoldChange) < (-logFC_threshold)), to_select := "yes"][order(log2FoldChange),]

    ## select the significant genes 
    significant_genes <- subset_data[to_select == "yes", Gene_Name]
  
    ## extract the name of the subset_data
    subdata_name <- names(DEgene_data_list[i])

    ## remove the underscore of the subsetdata_name to put in the upset plot
    subdata_name2 <- str_replace_all(
        string = subdata_name,
        pattern = "_vs_",
        replacement = " vs "
    )

    ## put the significant gene according with the subdata names
    significant_gene_list[[subdata_name2]] <- significant_genes
}

## take all the significant genes
significant_gene_vector <- unique(unlist(significant_gene_list, recursive = FALSE))

## take the significant genes in common amongst all the groups
significant_gene_vector_in_common <- Reduce(intersect, significant_gene_list)


##############################################
## Calculate the correlation between the genes of interest and the gene target
##############################################

## extract the expression data just associated with the significant gene vector
expression_data2 <- expression_data[ gene_id %in% significant_gene_vector, ]

## function for calculating the correlation
calculate_correlation <- function(
    gene_of_interest, # put the name of the gene of interest (gene name or gene symbols)
    expression_data_input
) {
    ## initilize the data table that will contains the correlation values
    correlation_dt_output <- data.table(
        "gene_of_interest" = character(),
        "target" = character(),
        "pearson" = numeric(),
        "pvalue" = numeric()
    )
    ## if the gene_of_interest is a ensembl id
    if (grepl(x = gene_of_interest, pattern = 'ENS') == FALSE) {
        gene_of_interest <- unlist(ensembl_id_gene_of_interest[gene_name == gene_of_interest, gene_id])
    }

    for (i_gene in seq(1, nrow(expression_data_input), 1)) {
        
        ## Extract the name of the gene target
        gene_target_name <- expression_data_input[i_gene, gene_id]

        ## extract the expression values of the gene on interest
        expression_gene_of_interest <- as.numeric(unlist(expression_data_input[grep(gene_id, pattern = gene_of_interest),]))

        ## extract the expression values of the gene target
        expression_gene_target <- as.numeric(unlist(expression_data_input[gene_id == gene_target_name,]))

        #expression_gene_of_interest <- 1:1000
        #expression_gene_target <- 1001:2000

        ## calculate the correlation bewteen the gene of interest and the gene target
        correlation_value <- cor(
            expression_gene_of_interest,
            expression_gene_target,
            method="pearson",
            use = "complete.obs"
        )

        ## do the statistical test
        pval <- cor.test(expression_gene_of_interest, expression_gene_target,  use = "complete.obs")$p.value
        
        ## create a new data table that will contain the results of iteration
        correlation_dt <- as.data.table(
            list(
                gene_of_interest,
                gene_target_name,
                correlation_value,
                pval
            )
        )[, setnames(
                .SD,
                c("V1", "V2", "V3", 'V4'),
                colnames(correlation_dt_output)
        )]

        #print(correlation_dt_output)
        ## concatenate the data together with the correlation_dt_output 
        correlation_dt_output <- rbind(
            correlation_dt_output,
            correlation_dt
        )
        print(correlation_dt_output)
    }

    return(correlation_dt_output)

}


## calculate the correlation for each genes
correlation_results <- foreach(i = gene_of_interest_vector, .combine='rbind', .packages='data.table') %dopar% calculate_correlation(i, expression_data = expression_data2)

#####
## calculate the pvalue adjusted by the bonferroni correction
#####

## extract the gene of interest
gene_of_interest_vector <- unlist(ensembl_id_gene_of_interest[, gene_id])

## copy the data table that contain the correlation results
correlation_data <- copy(correlation_results)

## calculate the pvalue adjusted by bonferroni correction
for (i_gene_of_interest in seq(1, length(gene_of_interest_vector), 1)) {

    ## extract the data asssociated with only one gene of interest
    subset_data <- correlation_data[gene_of_interest == gene_of_interest_vector[i_gene_of_interest],]

    ## based on the number of genes, correct the pvalues with bonferonni correction
    correlation_data[
        gene_of_interest == gene_of_interest_vector[i_gene_of_interest], 
        adj_pval :=  p.adjust(
            pvalue,
            n = nrow(subset_data),
            method = "bonferroni"
        )
    ]
}

##############################################
## write the results in a files
##############################################

## write the data table that contain
fwrite(
    correlation_data,
    snakemake@output[['correlation_results']],
    sep = ','
)
