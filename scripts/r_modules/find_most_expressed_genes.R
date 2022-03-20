
###################################################
###################################################


##### EXTRACT HE MOST EXPRESSED GENES FOR EACH GROUP AND ANALYSE THESE GENES


###################################################
###################################################


##### LOAD LIBRARIES
library(data.table)
library(stringr)
library(UpSetR)


###################################################
###################################################


##### LOAD THE DATA

## load the cluster groups data
cluster_group <- fread(
    file = snakemake@input[["cluster_group"]],
    sep = ","
)

## load the expression data
expression_data <- fread(
    file = snakemake@input[["expression_data"]],
    sep = ","
)

## load genes of interests
gene_id <- fread(
    file = snakemake@input[["gene_id"]],
    sep = ","
)

## load utils functions
source(snakemake@params[["utils_R"]])

## load the top gene thresholds
top_Expressed_Gene_threshold <- snakemake@params[["top_Expressed_Gene_threshold"]]

###################################################
###################################################


##### REFORMATE THE DATA


## formate the data ;extract the expression data related to the genes of interest ; transform the ensembl id to gene symbols and rename the gene name column by "Gene_Name"
expression_data <- expression_data[, rn := sapply(
    rn, function(x) str_split(
        x,
        pattern = "[.]",

    )[[1]][1]
)][
    rn %in% unlist(gene_id[, "gene_id"]),
][
    ,
    rn := ensemblID_to_geneSymbole(rn, gene_id)
][
    ,
    setnames(
        .SD,
        "rn",
        "Gene_Name"
    )
]


###################################################
###################################################


##### ANALYSIS


## function for extracting the most expressed gene for each groups
extract_top_genes <- function(
    CLUSTER_GROUP_INPUT, # data table of two columns ( "sample_id" and "cluster")
    EXPRESSION_DATA_INPUT, # expression data table : columns corresponding to id samples (+ one column called "Gen_Name corresponding to the gene names"), rows corresponding to genes,
    TOP_INPUT # corresponding to the number of genes that will be extracted
) {
    ## initialize the list of the top most expressed genes for each groups
    top_expressed_genes_list <- list()

    ## extract the cluster names from the cluster group data
    cluster_names <- unique(unlist(CLUSTER_GROUP_INPUT[, "cluster"]))

    ## for each cluster, extract the N first most expressed
    for (i_cluster_name in seq(1, length(cluster_names), 1)) {
        
        # print("----")
        # print(i_cluster_name)

        ## extract the sample ids associated with the cluster group
        cluster_sample <- unlist(CLUSTER_GROUP_INPUT[cluster == cluster_names[i_cluster_name], "sample_id", with = F])

        ## extract the expression data associated wutg the cluster sample
        cluster_expression_data <- EXPRESSION_DATA_INPUT[, append("Gene_Name", cluster_sample), with = F]

        ## do the mean for each genes
        cluster_expression_data <- cluster_expression_data[, expression_mean := rowMeans(.SD), .SDcols = cluster_sample][order(expression_mean, decreasing = T),]

        ## extract the gene name and the expression average
        cluster_gene_expression_mean <- (cluster_expression_data[, c("Gene_Name", "expression_mean"), with = F])

        ## extract the top genes
        top_genes <- list(unlist(cluster_gene_expression_mean[1:TOP_INPUT, "Gene_Name"]))

        ## Rename the list by cluster name
        names(top_genes) <- cluster_names[i_cluster_name]

        ## put into the list
        top_expressed_genes_list <- append(
            top_expressed_genes_list,
            top_genes,
        )
    }

    ## return the list
    return(top_expressed_genes_list)
}

## extract the top expressed genes
top_expressed_genes <- extract_top_genes(
    CLUSTER_GROUP_INPUT = cluster_group,
    EXPRESSION_DATA_INPUT = expression_data,
    TOP_INPUT = top_Expressed_Gene_threshold
)

## create the upset plot
upset_figure <- upset(
    fromList(top_expressed_genes),
    order.by = "freq",
    mainbar.y.label = "Number of Commun Genes",
    sets.x.label = "Number of Genes",
    
    point.size = 4,
    mb.ratio = c(0.55, 0.45),
    text.scale = 2
)

## save the uptset plot into the directory
svg(
    filename = snakemake@output[["upset_plot"]],
    width = 10,
    height = 7,
)
print(upset_figure)
dev.off()
