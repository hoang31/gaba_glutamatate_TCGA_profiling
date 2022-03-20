

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


## function for creating the immune figure
create_immune_stacked_barplot <- function(
    tool_name, # name of the name tool of interest
    directory_path
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
    colnames(subset_data) <- str_replace_all(colnames(subset_data), pattern = ' ', replacement = '_')
    
    ## retreive the correct name of each columns of interest
    colname_of_interest <- colnames(subset_data)[-(colnames(subset_data) %in% c('cluster'))]

    ## do the average for each cell type per groups
    average_data <- subset_data[, lapply(.SD, function(x) mean(x, na.rm = T)), by = c('cluster'), .SDcols = colname_of_interest]

    ## melt the data 
    average_data <- melt(
        average_data,
        id.vars = 'cluster',
        measure.vars = colname_of_interest
    )
    
    ## creathe the file path for the figure saving
    file_path <- paste(
        directory_path,
        "/",
        tool_name,
        '.png',
        sep = ''
    )
    if (tool_name == 'CIBERSORT') {

        ## generate the pie chart
        immune_plot <- ggplot(
            data = average_data,
            aes(
                title = tool_name,
                x = "",
                y = value,
                fill = variable
            )
        ) +
            labs(
                title = tool_name
            ) +
            scale_fill_manual(values = all_colors) +
            geom_bar(stat="identity") +
            coord_polar(theta = "y") +
            facet_grid(facets=. ~ cluster) + 
            theme_classic(base_size = 24) +
            theme(axis.text.x=element_blank())
    }
    else {

        ## generate the stacked barchart
        immune_plot <- ggplot(
            data = average_data,
            aes(
                title = tool_name,
                x = cluster,
                y = value,
                fill = variable
            )
        ) +
            labs(
                title = tool_name,
                y = "Absolute Quantification"
            ) +
            scale_fill_manual(values = all_colors) +
            geom_bar(stat="identity") +
            #coord_polar(theta = "y") +
            #facet_grid(facets=. ~ cluster) + 
            theme_classic(base_size = 24)
    }

    ggsave(
        file_path,
        plot = immune_plot,
        device = 'png',
        width = 18,
        height = 8
    )

    return(immune_plot)
}


########################################
## Write all the immune pltos into files
########################################

## create the directory
dir.create(snakemake@output[["immune_composition_figures_dir"]])

## set the data of which tools we will used
tools <- c(
    'CIBERSORT',
    'CIBERSORT-ABS',
    'QUANTISEQ',
    'EPIC',
    'TIMER',
    'XCELL'
)

## color to use in the plot
all_colors <- c(
    brewer.pal(9, "Set1"),
    brewer.pal(8, "Dark2"),
    brewer.pal(8, "Set2"),
    brewer.pal(12, "Paired"),
    brewer.pal(9, "Set1")
)


## create the immune plots for every selected tools
immune_plot_list <- lapply(
    tools,
    function(x) create_immune_stacked_barplot(x, directory_path = snakemake@output[["immune_composition_figures_dir"]])
)

## create a combined plot composed of the plots generated from cibersort and cibersort-abs
combined_plot <- ggarrange(immune_plot_list[[1]], immune_plot_list[[2]], ncol = 1)
ggsave(
    paste(
        snakemake@output[["immune_composition_figures_dir"]],
        "/CIBERSORT_combined_plot.svg",
        sep = ""
    ),
    plot = combined_plot,
    device = 'svg',
    width = 20,
    height = 18
)


########################################
## Do the statistical analysis between each clusters
########################################

#####
## function extracting the data related to the tool of interest
#####

extract_immune_data <- function(
    tool_name, # name of the name tool of interest
    immune_data
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
    colname_of_interest <- grep(x = colnames(immune_data), pattern = tools_regex, value = T)

    ## extract the data associated with the column name of interest
    subset_data <- immune_data[,c('sample_id', 'cluster', colname_of_interest), with = F]

    ## transform all the values to numeric
    subset_data[, (colname_of_interest) := lapply(.SD, as.numeric), .SDcols = colname_of_interest]

    ## rename the column name
    colnames(subset_data) <- str_replace_all(colnames(subset_data), pattern = tools_regex, replacement = '')
    colnames(subset_data) <- str_replace_all(colnames(subset_data), pattern = ' ', replacement = '_')
    
    ## retreive the correct name of each columns of interest
    colname_of_interest <- colnames(subset_data)[-(colnames(subset_data) %in% c('cluster'))]

    ## return the data related to the tool of interest
    return(subset_data)
}

#####
## function for the stastistical test for each immune cell type
#####

statistical_analysis <- function(
    immune_data # immune data related to one tool of interest (output of the extract_immune_data)
) {
    
    ## transform the immune data to data table
    immune_data <- as.data.table(immune_data)

    ## initialize the data table that will contain the statistical results
    statistical_results_dt <- data.table(
        "cell_type" = character(),
        "cluster1" = character(),
        "cluster2" = character(),
        "pvalue" = numeric()
    )

    ## identification of the cluster names
    cluster_name_vector <- unique(unlist(immune_data[, cluster]))

    ## do the paired combination of the cluster name vector
    cluster_combination <- combn(
        x = cluster_name_vector,
        m = 2
    )
    
    ## identification of the cell types contained by the immune data
    cell_type_vector <- colnames(immune_data)[!(colnames(immune_data) %in% c("sample_id", "cluster"))]

    ## for each cell type
    for (i_cell_type in seq(1, length(cell_type_vector), 1)) {

        ## retreive the cell type 
        cell_type <- cell_type_vector[i_cell_type]
        
        for (i_combination in seq(1, ncol(cluster_combination), 1)) {
            
            ## extract the clusters of interest
            cluster1 <- cluster_combination[1, i_combination]
            cluster2 <- cluster_combination[2, i_combination]

            ## retreive the values related to each cluster and cell type of interest
            cluster1_values <- unlist(immune_data[cluster == cluster1, cell_type, with = F])
            cluster2_values <- unlist(immune_data[cluster == cluster2, cell_type, with = F])

            ## do the statistical test
            pval <- wilcox.test(
                cluster1_values,
                cluster2_values
            )$p.value

            ## put the results into a data table
            results_dt <- data.table(
                "cell_type" = cell_type,
                "cluster1" = cluster1,
                "cluster2" = cluster2,
                "pvalue" = pval
            )

            ## put the results into the statistical_results_dt
            statistical_results_dt <- rbind(
                statistical_results_dt,
                results_dt
            )

        }
        
    }

    ## return the statistical results
    return(statistical_results_dt)

}


#####
## call the functions
#####

## set the data of which tools we will used
 tools <- c(
    'CIBERSORT',
    'CIBERSORT-ABS',
    'QUANTISEQ',
    'EPIC',
    'TIMER',
    'XCELL'
 )

## extract all the data related to the tools of interest into a list
immune_data_list <- lapply(
    tools,
    function(x) extract_immune_data(
        tool_name = x,
        immune_data = immune_composition_data
    )
)

## do the statistical analysis for all the immune data contained in the immune_data_list
statistical_results_list <- lapply(
    immune_data_list,
    function(x) statistical_analysis(immune_data = x)
)

## rename the list
names(statistical_results_list) <- tools

## write the statistical data
for (i_stat_results in seq(1, length(statistical_results_list), 1)) {

    ## generate the path of the file
    file_path <- paste(
            snakemake@output[["immune_composition_figures_dir"]],
            "/",
            names(statistical_results_list[i_stat_results]),
            "_pval",
            ".csv",
            sep = ""
    )
    
    ## write the data table that contain the stats into a file
    fwrite(
        x = as.data.table(statistical_results_list[i_stat_results]),
        file = file_path,
        sep = ","
    )
}



##########
## Statistical analysis between cell type quantification (cibersort abs)
##########


## extract data related to the cibersort abs 
cibersort_abs_data <- extract_immune_data(
    'CIBERSORT-ABS',
    immune_composition_data

)

## extract the cell type name contained in the data
cell_type_vector <- colnames(cibersort_abs_data)[!(colnames(cibersort_abs_data) %in% c('sample_id', 'cluster'))]

## melt the data
cibersort_abs_data <- melt(
    cibersort_abs_data,
    id.vars = c('sample_id', 'cluster'),
    measure.vars = cell_type_vector
)

## do the mean by sample
cibersort_abs_average <- cibersort_abs_data[, .(average = mean(value)), by = c("sample_id", "cluster")]

## do the statistical test
stat_res <- statistical_analysis(
    cibersort_abs_average
)

## do the bferonni correction
stat_res <- stat_res[, adj_pvalue := p.adjust(pvalue, method = "bonferroni", n = nrow(stat_res))][order(adj_pvalue, decreasing = F),]


## write the statistical test
fwrite(
    stat_res,
    paste(
            snakemake@output[["immune_composition_figures_dir"]],
            "/",
            "pval_cibersortabs_quantification_per_cluster",
            ".csv",
            sep = ""
    ),
    sep = ","
)