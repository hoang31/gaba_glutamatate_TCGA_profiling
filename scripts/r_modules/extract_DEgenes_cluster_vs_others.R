
###############################################
## Extract the DE genes from the cluster vs others analysis
###############################################


###############################################
## Load the libraries
###############################################

library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(UpSetR)
library(EnhancedVolcano)
library(egg)


###############################################
## Load the data 
###############################################

## load the cluster data
clusters <- fread(
    snakemake@input[["cluster_group"]],
    sep = ","
)

## retreive the cluster names of the data
cluster_names <- unique(unlist(clusters[, cluster]))

## load the genes related to GABA glutamate and calcium pathway
gene_vector <- unique(
    fread(
        snakemake@input[["genes_highExpressed"]],
        sep = ","
    )[, gene_name]
)

## Set the pvalue cutoff and the fold change cutoff
pval_cutoff <- 0.001
FCcutoff <- 1

## create the directory that will contain the figures
dir.create(snakemake@output[["DEgene_figures"]])


###############################################
## Analysis
###############################################


##########
## Load the DE results
##########

## load the directory
DEgenes_directory <- snakemake@input[["DEgene_results_dir"]]

## retreive the file names
file_names <- list.files(DEgenes_directory)

## retreive only the file names that describe just the filtered results
file_names <- grep(file_names, pattern = "results", value = TRUE)

## create the file path from the file names
file_path <- paste(
    snakemake@input[["DEgene_results_dir"]],
    file_names,
    sep = "/"
)

## read all the files using the path files
merged_file <- lapply(
    X = file_path,
    function(x) {
        res <- fread(
            input = x,
            sep = ","
        )
        return(res)
    }
)

## transform the names for each file
file_names <- str_replace_all(file_names, pattern = "results_", replacement = "")
file_names <- str_replace_all(file_names, pattern = ".csv", replacement = "")

## rename the names of each data removing the filtered_ word from the name
names(merged_file) <- file_names


##########
## Load the KEGG pathway genes
##########

## load the KEGG pathways genes paths
KEGG_genes_file_path <- list.files(
    snakemake@input[["kegg_pathway_genes_directory"]],
    full.names = T
)

## load the genes
KEGG_genes_list <- lapply(
    X = KEGG_genes_file_path,
    function(x) {
        list_output <- scan(x, what = character(), sep = "\n")
        return(list_output)
    }
)

## take the KEGG pathway names
KEGG_genes_file_names <- list.files(
    snakemake@input[["kegg_pathway_genes_directory"]],
    full.names = F
)

## remove the .txt extension
KEGG_genes_file_names <- str_replace(
    KEGG_genes_file_names,
    pattern = ".txt",
    replacement = ""
)

## remane the kegg gene list with the kegg pathway names
names(KEGG_genes_list) <- KEGG_genes_file_names

## initialize a new list
KEGG_genes_list2 <- list()

## merge the common genes into a same pathways
KEGG_genes_list2[["calcium"]] <- unique(c(
    KEGG_genes_list[["CALCIUM_endocrine"]],
    KEGG_genes_list[["CALCIUM_signaling"]]
))

KEGG_genes_list2[["glutamate_pathway"]] <- unique(c(
    KEGG_genes_list[["GLUTAMATE_synapse"]],
    KEGG_genes_list[["GLUTAMATE_metabolism"]]
))

KEGG_genes_list2[["gaba_pathway"]] <- unique(c(
    KEGG_genes_list[["GABA_synapse"]]
))

KEGG_genes_list2[["gaba_gluta_pathway"]] <- intersect(
    KEGG_genes_list2[["glutamate_pathway"]],
    KEGG_genes_list2[["gaba_pathway"]]
)

##########
## GENERATE THE VOLCANOPLOTS FOR EACH RESULTS OF DESEQ2
##########


## for each data in the merged_file
for (i in seq(1, length(merged_file), 1)) {


    ## extract the subset_data of the merged_file list
    subset_data <- as.data.table(merged_file[[i]])

    ## extract only the results related to the gene of interest
    subset_data <- subset_data[Gene_Name %in% gene_vector,]

    ## extract the name of the subset_data
    subdata_name <- names(merged_file[i])

    ## transform to data frame and rename the row names with the gene names
    setDF(subset_data, rownames = unlist(subset_data[,"Gene_Name"]))

    ## Initialize the table data for the color associated with each kegg pathway and genes
    dt_kegg_pathway_colors <- data.table(
        "Gene_Name" = rownames(subset_data)
    )

    ## color palette
    colors <- brewer.pal(9, "Set1")

    ## add the color in the data table for each pathway
    for (i in seq(1, length(KEGG_genes_list2), 1)) {
        dt_kegg_pathway_colors[Gene_Name %in% unlist(KEGG_genes_list2[[i]]), "pathway" := names(KEGG_genes_list2[i])]
        dt_kegg_pathway_colors[Gene_Name %in% unlist(KEGG_genes_list2[[i]]), "color" := colors[i]]
    }

    ## transform the subset data into data table for subsetting easier
    subset_data2 <- copy(subset_data)
    setDT(subset_data2)

    ## extract genes which are below the pvalue and fc cutoff
    ns_genes <- unlist(subset_data2[((log2FoldChange >= -FCcutoff) & (log2FoldChange <= FCcutoff)) | padj > pval_cutoff, "Gene_Name"])
    
    ## color the genes which are not significant in grey and put "NS" as label
    dt_kegg_pathway_colors[Gene_Name %in% ns_genes, "pathway" := "NS"]
    dt_kegg_pathway_colors[Gene_Name %in% ns_genes, "color" := "grey"]

    ## create the color vectors that will be used in the volcano plot generation
    KEGG_color <- dt_kegg_pathway_colors[["color"]]
    names(KEGG_color) <- dt_kegg_pathway_colors[["pathway"]]

    ## generate the title of the plot
    title <- str_replace_all(
        subdata_name,
        pattern = "_",
        replacement = " "
    )

    ## generate the Volcanoplot
    volca_plot <- EnhancedVolcano(
        toptable = subset_data,
        lab = rownames(subset_data),
        selectLab = c(''),
        x = "log2FoldChange",
        y = "padj",

        title = "",
        pCutoff = pval_cutoff,
        FCcutoff = FCcutoff,

        xlim = c(-5, 5),
        colAlpha = 0.99,
        pointSize = 6,
        labSize = 5,

        colCustom = KEGG_color
    )


    ## add legend
    volca_plot <- volca_plot + 
    theme_minimal(base_size = 24) +
    labs(color='Pathways') +
    labs(subtitle = title)


    ## generate the file path for saving
    file_path <- paste(
        snakemake@output[["DEgene_figures"]],
        "/volcano_plot_",
        subdata_name,
        ".svg",
        sep = ""
    )

    ## save the volcano plot in the diretory
    svg(
        filename = file_path,
        width = 14,
        height = 8
    )
    #png(
    #    filename = file_path,
    #    width = 1000,
    #    height = 600
    #)
    print(volca_plot)
    dev.off()
    

}

##########
## GENERATE UPSETPLOTS
##########

## initialize the list which contains the significant pvalues for each subset of the merged_data
significant_gene_list <- list()

## initialize the p value and the logfoldchange thresolds
pval_threshold <- pval_cutoff
logFC_threshold <- FCcutoff

## for each subset of the merged_data, extract the genes associated with a pvalue lower than pval_threshold and put them into the significant_gene_list
for (i in seq(1, length(merged_file), 1)) {

    # print("###############################")

    ## extract the subset_data of the merged_file list
    subset_data <- as.data.table(merged_file[[i]])

    ## extract only the results related to the gene of interest
    subset_data <- subset_data[Gene_Name %in% gene_vector,]

    ## we have a problem when we are selecting the genes depending of the log fold change so I foound a workaround
    subset_data[, log2FoldChange := as.numeric(log2FoldChange)]
    subset_data <- (subset_data[(padj < pval_threshold),])
    subset_data <- subset_data[((log2FoldChange) > logFC_threshold), to_select := "yes"][order(log2FoldChange),]
    subset_data <- subset_data[((log2FoldChange) < (-logFC_threshold)), to_select := "yes"][order(log2FoldChange),]

    ## select the significant genes 
    significant_genes <- subset_data[to_select == "yes", Gene_Name]
  
    ## extract the name of the subset_data
    subdata_name <- names(merged_file[i])

    ## remove the underscore of the subsetdata_name to put in the upset plot
    subdata_name2 <- str_replace_all(
        string = subdata_name,
        pattern = "_vs_",
        replacement = " vs "
    )

    ## put the significant gene according with the subdata names
    significant_gene_list[[subdata_name2]] <- significant_genes
}

## generate the upset plot with the significant gene list
upset_plot <- upset(
    fromList(significant_gene_list),
    order.by = "freq",
    nsets = length(significant_gene_list),
    mainbar.y.label = "Number of Common Genes",
    sets.x.label = "Number of Genes",
    
    point.size = 4,
    mb.ratio = c(0.55, 0.45),
    text.scale = 2
)


## generate the file path for saving
file_path <- paste(
    snakemake@output[["DEgene_figures"]],
    "/upset_plot_DEgenes",
    ".svg",
    sep = ""
)

## save the uptset plot into the directory
svg(
    filename = file_path,
    width = 16,
    height = 7,
)
print(upset_plot)
dev.off()
