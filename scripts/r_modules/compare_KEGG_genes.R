
###################################################
###################################################


##### EXTRACT THE MOST EXPRESSED GENES FOR EACH GROUP AND ANALYSE THESE GENES


###################################################
###################################################


##### LOAD LIBRARIES
library(data.table)
library(stringr)
library(UpSetR)


###################################################
###################################################


##### LOAD THE DATA

## load the directory containing the kegg pathways genes
kegg_pathway_genes_dir <- snakemake@input[["kegg_pathway_genes_directory"]]


## load the directory files
kegg_pathway_genes_files <- list.files(kegg_pathway_genes_dir)

## add the path
kegg_pathway_genes_files <- sapply(kegg_pathway_genes_files, function(x) paste(
    kegg_pathway_genes_dir,
    x,
    sep = "/"
))

## put all to the same list
kegg_pathway_genes_files_list <- lapply(
    kegg_pathway_genes_files,
    function(x) unlist(fread(
        x,
        sep = ","
    ))
)

## change the names removing the ".txt" extension
names(kegg_pathway_genes_files_list) <- sapply(
    names(kegg_pathway_genes_files_list),
    function(x) str_split(
        x,
        pattern = "[.]"
    )[[1]][1]
)


## count the number of unique genes 
length(unique(unlist(kegg_pathway_genes_files_list, recursive = T)))


###################################################
###################################################


###### GENERATE THE UPSET PLOT

## create the upset plot
upset_figure <- upset(
    fromList(kegg_pathway_genes_files_list),
    order.by = "freq",
    mainbar.y.label = "Number of Commun Genes",
    sets.x.label = "Number of Genes",
    show.numbers = 'yes',
    point.size = 4,
    mb.ratio = c(0.55, 0.45),
    text.scale = 2
)
## save the uptset plot into the directory
svg(
    filename = snakemake@output[["upset_plot_genes_of_interest"]],
    width = 10,
    height = 7,
)
print(upset_figure)
dev.off()
