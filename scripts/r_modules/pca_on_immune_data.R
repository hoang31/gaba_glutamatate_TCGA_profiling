

########################################
## Infer immune composition for each sample
########################################


########################################
## Load the libraries
########################################


library(data.table)
library(tidyverse)
library(RColorBrewer)
library(egg)
library(factoextra)
library(FactoMineR)


########################################
## load the data
########################################

## load the immune compisition data 
immune_composition_data <- fread(
    snakemake@input[['immune_composition_data']],
    sep = ','
)

## load the clinical data
clinical_data <- fread(
    snakemake@input[['clinical_data']],
    sep = ','
)

## load the cluster data
cluster_group <- fread(
    snakemake@input[['cluster_group']],
    sep = ','
)


########################################
## Analysis
########################################


##########
## Extract the immune data of the sample of interest
##########

## extract the sample id that we are intereted to
sample_id_vector <- unlist(clinical_data[, Sample_ID])

## reformate the data removing the letter at the end
sample_id_vector <- sapply(sample_id_vector, function(x) str_replace_all(x, pattern = 'A$', replacement = ''))
sample_id_vector <- sapply(sample_id_vector, function(x) str_replace_all(x, pattern = 'B$', replacement = ''))
sample_id_vector <- sapply(sample_id_vector, function(x) str_replace_all(x, pattern = 'C$', replacement = ''))

## extract the data of interest
immune_composition_data <- immune_composition_data[cell_type %in% sample_id_vector,]

## formate the sample id
immune_composition_data[, cell_type := sapply(cell_type, function(x) paste(str_split(x, pattern = '[-]')[[1]][1:3], collapse = '-'))]

## merge the cluster and the immune data together
immune_composition_data <- merge(
    cluster_group,
    immune_composition_data,
    by.x = 'sample_id',
    by.y = 'cell_type',
    sort = F,
    all.x = T
)

##########
## Function
##########

## function for extracting the data associated with a tool
extract_immune_data <- function(
    tool_name # name of the name tool of interest
) {

    ## create the regex from the tools names
    tools_regex <- paste(
        "_",
        tool_name,
        '$',
        sep = '',
        collapse = '|'
    )

    ## extract the column names of interest
    colname_of_interest <- grep(x = colnames(immune_composition_data), pattern = tools_regex, value = T)

    ## extract the data associated with the column name of interest
    subset_data <- immune_composition_data[,c('sample_id', 'cluster', colname_of_interest), with = F]
    subset_data <- immune_composition_data[,c('cluster', colname_of_interest), with = F]

    ## transform all the values to numeric
    subset_data[, (colname_of_interest) := lapply(.SD, as.numeric), .SDcols = colname_of_interest]

    ## rename the column name
    colnames(subset_data) <- str_replace_all(colnames(subset_data), pattern = tools_regex, replacement = '')
    #colnames(subset_data) <- str_replace_all(colnames(subset_data), pattern = ' ', replacement = '_')
    
    ## retreive the correct name of each columns of interest
    colname_of_interest <- colnames(subset_data)[-(colnames(subset_data) %in% c('cluster'))]

    ## return the data
    return(subset_data)
}


## function for performing the PCA and generate related figures
perform_pca <- function(
    tool_name # string that corresponds to the tool names
) {

    ## extract the data related to the tool and remove the na values
    tool_data <- na.omit(
            extract_immune_data(
            tool_name = tool_name
        )
    )

    #tool_data <- tool_data[cluster %in% c("IDHmut_CODEL", "IDHwt"), ]

    ## perform the pca
    res.pca <- prcomp(
        tool_data[, !("cluster"), with = F],
        scale = FALSE
    )

    ## generate the barchart that show the contribution for each coponent
    contribution_barchart <- fviz_eig(res.pca)

    ## generate the plot related to the sample data
    sample_plot <- fviz_pca_ind(
        res.pca,
        col.ind = "cos2", # Colorer par le cos2
        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
        repel = TRUE     
    ) + 
        theme_classic(base_size = 24)

    ## generate the file path of the figure to save
    file_path <- paste(
        snakemake@output[["pca_figure_dir"]],
        "/",
        tool_name,
        "_",
        "sample_plot",
        ".svg",
        sep = ""
    )

    ## save the plot of interest
    ggsave(
        sample_plot,
        filename = file_path,
        device = "svg",
        height = 8,
        width = 12
    )

    ## generate variable plot
    variable_plot <- fviz_pca_var(
        res.pca,
        col.var = "contrib", 
        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
        repel = TRUE,
        labelsize = 7,
        arrowsize = 1
    ) + 
        theme_classic(base_size = 24)

    ## generate the file path of the figure to save
    file_path <- paste(
        snakemake@output[["pca_figure_dir"]],
        "/",
        tool_name,
        "_",
        "variable_plot",
        ".svg",
        sep = ""
    )

    ## save the plot of interest
    ggsave(
        variable_plot,
        filename = file_path,
        device = "svg",
        height = 10,
        width = 12
    )

    #fviz_pca_biplot(
    #    res.pca,
    #    repel = TRUE,
    #    col.var = "#2E9FDF", 
    #    col.ind = "#696969"  
    #)


    ## retrieve the variable contribution for each component
    contribution <- PCA(
        tool_data[, !("cluster"), with = F],
        scale.unit = FALSE
    )

    contribution <- get_pca_var(contribution)

    print("----- variable contribution -----")
    print(contribution$contrib)

 
    ## contribution for each variable for the first two component
    contribution_1axe <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10) +
        theme_classic(base_size = 24) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    contribution_2axe <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10) +
        theme_classic(base_size = 24) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ## generate the file path of the figure to save
    file_path <- paste(
        snakemake@output[["pca_figure_dir"]],
        "/",
        tool_name,
        "_",
        "contribution_axe1",
        ".svg",
        sep = ""
    )
    file_path2 <- paste(
        snakemake@output[["pca_figure_dir"]],
        "/",
        tool_name,
        "_",
        "contribution_axe2",
        ".svg",
        sep = ""
    )

    ## save the plot of interest
    ggsave(
        contribution_1axe,
        filename = file_path,
        device = "svg",
        height = 10,
        width = 12
    )
    ggsave(
        contribution_2axe,
        filename = file_path2,
        device = "svg",
        height = 10,
        width = 12
    )

    #####
    ## generate the variable contribution plot
    #####

    sample2_plot <- fviz_pca_ind(
        res.pca,
        label = 'none',
        col.ind = as.factor(unlist(tool_data[, c("cluster"), with = F])), 
        palette = c("#00AFBB",  "#FC4E07", "#fcf807", "#3dac21"),
        #palette = c("#ffffff00",  "#ffffff00", "#fcf807", "#3dac21"),

        addEllipses = TRUE,
        #ellipse.type = "confidence",
        legend.title = "Groups",
        repel = TRUE
    ) + 
        theme_classic(base_size = 24)

    ## generate the file path of the figure to save
    file_path <- paste(
        snakemake@output[["pca_figure_dir"]],
        "/",
        tool_name,
        "_",
        "sample2_plot",
        ".svg",
        sep = ""
    )

    ## save the plot of interest
    ggsave(
        sample2_plot,
        filename = file_path,
        device = "svg",
        height = 8,
        width = 12
    )

}

########################################
## Save the pca figures into the directory
########################################

## create the directory that will contain the figures
dir.create(snakemake@output[["pca_figure_dir"]])

## generate the pca on the CIBERSORTx tools
perform_pca(
    tool_name = 'CIBERSORT-ABS'
)


## generate the pca on the CIBERSORT tool
perform_pca(
    tool_name = 'CIBERSORT'
)


