
########################################
## Do the enrichment analysis of the genes that correlate with the gene of interest
########################################

########################################
## Load the libraries
########################################

library(data.table)
library(tidyverse)
library(clusterProfiler)
library(GO.db)

options(width=system("tput cols", intern=TRUE))

########################################
## Load the data
########################################

## load the correlation data
correlation_data <- fread(
    snakemake@input[['correlation_results']],
    sep = ","
)

## initialize the directory that will contain the output
dir.create(snakemake@output[['dir_correlated_gene_enrichment_results']])

## set the cutoff for considering the significativity
pval_cutoff <- 0.001

## load external scripts
source(snakemake@params[["utils_R"]])

########################################
## Preparation of the data
########################################

## add into the correlation data some information about the correlation
correlation_data[pearson < 0, correlation := 'negative']
correlation_data[pearson > 0, correlation := 'positive']
correlation_data[adj_pval > pval_cutoff, correlation := 'ns']

## display the number of genes that correlate negatively or negatively
print(table(correlation_data[, correlation]))

## put the same gene name for every genes
correlation_data[, gene_of_interest := "gene1"]

## extract the genes of interest
genes_of_interest_vector <- unique(unlist(correlation_data[,gene_of_interest]))

## initialize the list that will contain all the enrichment results
GO_enrichment_result_list <- list()

## for each gene of interest, extract all the genes that correlate with
for (i_gene_of_interest in seq(1, length(genes_of_interest_vector), 1)) {

    ## extract the data of the gene of interest
    correlation_subset_data <- correlation_data[gene_of_interest == genes_of_interest_vector[i_gene_of_interest],]

    ## print the informations about the distribution of the correlated genes
    print(table(correlation_subset_data[,correlation]))
    print(table(correlation_subset_data[,correlation])/nrow(correlation_subset_data))

    ## extract only negative or positive correlated genes
    positive_correlated_gene_vector <- unique(unlist(correlation_subset_data[correlation == "positive",target]))
    negative_correlated_gene_vector <- unique(unlist(correlation_subset_data[correlation == "negative",target]))
    ns_correlated_gene_vector <- unique(unlist(correlation_subset_data[correlation == "ns",target]))

    ## remove the part after the point of the ensembl id
    positive_correlated_gene_vector <- sapply(
        positive_correlated_gene_vector, 
        function(x) str_split(x, pattern = '[.]')[[1]][1]
    )
    negative_correlated_gene_vector <- sapply(
        negative_correlated_gene_vector, 
        function(x) str_split(x, pattern = '[.]')[[1]][1]
    )
    ns_correlated_gene_vector <- sapply(
        ns_correlated_gene_vector, 
        function(x) str_split(x, pattern = '[.]')[[1]][1]
    )

    ## transform the gene of interest ensembl id to symbol
    #gene_of_interest_symbol <- bitr(
    #    genes_of_interest_vector[i_gene_of_interest],
    #    fromType = "ENSEMBL",
    #    toType = c("SYMBOL", "ENTREZID"),
    #    OrgDb = 'org.Hs.eg.db'
    #)[, 'SYMBOL']
    
    gene_of_interest_symbol <- "gene1"


    ## transform the ensembl id to other id types
    positive_correlated_gene_vector <- bitr(
        positive_correlated_gene_vector,
        fromType = "ENSEMBL",
        toType = c("SYMBOL", "ENTREZID"),
        OrgDb = 'org.Hs.eg.db'
    )[, 'ENTREZID']
    
    negative_correlated_gene_vector <- bitr(
        negative_correlated_gene_vector,
        fromType = "ENSEMBL",
        toType = c("SYMBOL", "ENTREZID"),
        OrgDb = 'org.Hs.eg.db'
    )[, 'ENTREZID']


    ns_correlated_gene_vector <- bitr(
        ns_correlated_gene_vector,
        fromType = "ENSEMBL",
        toType = c("SYMBOL", "ENTREZID"),
        OrgDb = 'org.Hs.eg.db'
    )[, 'ENTREZID']

    
    ##########
    ## do the KEGG enrichment analysis on the positive and the negative correlated genes
    ##########

    kegg_res_positive <- as.data.table(summary(enrichKEGG(
        gene = positive_correlated_gene_vector,
        organism = 'hsa',
        pvalueCutoff = 0.05

    )))[order(Count, decreasing = T),]

    kegg_res_negative <- as.data.table(summary(enrichKEGG(
        gene = negative_correlated_gene_vector,
        organism = 'hsa',
        pvalueCutoff = 0.05
    )))[order(Count, decreasing = T),]

    kegg_res_ns <- as.data.table(summary(enrichKEGG(
        gene = ns_correlated_gene_vector,
        organism = 'hsa',
        pvalueCutoff = 0.05
    )))[order(Count, decreasing = T),]

    ## write the kegg results into files
    fwrite(
        kegg_res_positive,
        paste(
            snakemake@output[['dir_correlated_gene_enrichment_results']],
            '/kegg_positive.csv',
            sep = ''

        ),
        sep = ','
    )

    fwrite(
        kegg_res_negative,
        paste(
            snakemake@output[['dir_correlated_gene_enrichment_results']],
            '/kegg_negative.csv',
            sep = ''
        ),
        sep = ','
    )

    fwrite(
        kegg_res_ns,
        paste(
            snakemake@output[['dir_correlated_gene_enrichment_results']],
            '/kegg_ns.csv',
            sep = ''
        ),
        sep = ','
    )

    ##########
    ## do the GO enrichment analysis on the positive and the negative correlated genes
    ##########

    GO_terms <- c(
        'CC',
        'MF',
        'BP'
    )

    for (i_GO_terms in seq(1, length(GO_terms), 1)) {

        print(GO_terms[i_GO_terms])
        go_res_positive <- enrichGO(
            gene = positive_correlated_gene_vector,
            OrgDb = 'org.Hs.eg.db',
            ont = GO_terms[i_GO_terms],
            pAdjustMethod = "bonferroni",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.05,
            readable = TRUE
        )

        go_res_negative <- enrichGO(
            gene = negative_correlated_gene_vector,
            OrgDb = 'org.Hs.eg.db',
            ont = GO_terms[i_GO_terms],
            pAdjustMethod = "bonferroni",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.05,
            readable = TRUE
        )

        go_res_ns <- enrichGO(
            gene = ns_correlated_gene_vector,
            OrgDb = 'org.Hs.eg.db',
            ont = GO_terms[i_GO_terms],
            pAdjustMethod = "bonferroni",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.05,
            readable = TRUE
        )

        ## formate the results
        go_res_positive <- as.data.table(summary(go_res_positive))[order(Count, decreasing = T),][, data_type := 'positive']
        go_res_negative <- as.data.table(summary(go_res_negative))[order(Count, decreasing = T),][, data_type := 'negative']
        go_res_ns <- as.data.table(summary(go_res_ns))[order(Count, decreasing = T),][, data_type := 'ns']

        ## create the paths for saving the files
        file_path_positive <- paste(
            snakemake@output[['dir_correlated_gene_enrichment_results']],
            '/',
            gene_of_interest_symbol,
            '_GO_',
            GO_terms[i_GO_terms],
            '_positive_correlation.csv',
            sep = ""
        )
        file_path_negative <- paste(
            snakemake@output[['dir_correlated_gene_enrichment_results']],
            '/',
            gene_of_interest_symbol,
            '_GO_',
            GO_terms[i_GO_terms],
            '_negative_correlation.csv',
            sep = ""
        )
        file_path_ns <- paste(
            snakemake@output[['dir_correlated_gene_enrichment_results']],
            '/',
            gene_of_interest_symbol,
            '_GO_',
            GO_terms[i_GO_terms],
            '_ns_correlation.csv',
            sep = ""
        )

        ## write the data
        fwrite(
            go_res_positive,
            file_path_positive,
            sep = ','
        )
        fwrite(
            go_res_negative,
            file_path_negative,
            sep = ','
        )
        fwrite(
            go_res_ns,
            file_path_ns,
            sep = ','
        )

        ## add the GO enrichment results in the output list
        GO_enrichment_result_list[[GO_terms[i_GO_terms]]]['positive'] <- go_res_positive
        GO_enrichment_result_list[[GO_terms[i_GO_terms]]]['negative'] <- go_res_negative
        GO_enrichment_result_list[[GO_terms[i_GO_terms]]]['ns'] <- go_res_ns

    }
}


########################################
## generate the figures
########################################

## copy the correlation data for modification
correlation_data2 <- copy(correlation_data)
correlation_data2[adj_pval > pval_cutoff, correlation := 'ns']

## order the levels for the ledend
correlation_data2$correlation <- factor(
    correlation_data2$correlation,
    levels = c(
        'positive',
        'negative',
        'ns'
    ),
)

#####
## Plot that described the pearson correlation depending to the p values
#####

## create the plot of the significant correlated genes 
correlated_gene_plot <- ggplot(
    correlation_data2,
    aes(
        x = pearson,
        y = -log10(adj_pval)
    )
) +
    geom_point(aes(color = correlation), size = 3) +
    #geom_point(shape = 1, size = 7, colour = "black") +
    scale_color_manual(
        values = c(
            '#ec8585',
            '#6262d6',
            'grey'
        )
    ) +
    #scale_x_continuous(
        #breaks = seq(-1, 1, by = 0.1)
    #) +
    theme_classic(
        base_size = 24
    ) +
    labs(
        x = 'Pearson Coefficient',
        y = '-log10(adjusted pvalue)'
    )

ggsave(
    correlated_gene_plot,
    filename = paste(
        snakemake@output[['dir_correlated_gene_enrichment_results']],
        '/correlated_gene_plot.png',
        sep = ''
    ),
    device = 'png',
    height = 8,
    width = 16
)

#####
## Generate the density plot of the pearson coefs
#####

## create the density plot of the correlation
correlation_density_plot <- ggplot(
    data = correlation_data2,
    aes(x = pearson)
) +
    geom_density(size = 3) +
    theme_classic(
        base_size = 24
    ) 

ggsave(
    correlation_density_plot,
    filename = paste(
        snakemake@output[['dir_correlated_gene_enrichment_results']],
        '/correlated_gene_density_plot.png',
        sep = ''
    ),
    device = 'png',
    height = 8,
    width = 16
) 
