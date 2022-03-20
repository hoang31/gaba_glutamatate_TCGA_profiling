
#############################################
## Generate the figures of correlation analysis between methylation and expression data
#############################################

#############################################
## Load the libraries
#############################################


library(data.table)
library(tidyverse)
library(RColorBrewer)


## set the color palette
color_palette <- brewer.pal(9, "Set1")


#############################################
## Load the data
#############################################


## load the utils module
source(snakemake@params[['utils']])

## read the data table that contain the methylation information
methylation_information_dt <- fread(
    snakemake@input[["methylation_information"]],
    sep = ","
)

## load the methylation data
methylation_data_path_vector <- list.files(
    snakemake@input[["DMgenes_analysis_dir"]],
    full.names = T
)
methylation_data_dt <- fread(
    grep(methylation_data_path_vector, pattern = "methylation_data_filtered.csv", value = T),
    sep = ","
)

## remove healthy sample from the methylation data
healthy_sample_vector <- c(
    "TCGA-06-AABW",
    "TCGA-74-6573"
)
methylation_data_dt <- methylation_data_dt[!(sample_id %in% healthy_sample_vector),]

## load the expression data
expression_data_dt <- fread(
    snakemake@input[["expression_data"]],
    sep = ","
)

## load the genes related to GABA glutamate and calcium pathway
neurotransmission_gene_vector <- unique(
    fread(
        snakemake@input[["genes_highExpressed"]],
        sep = ","
    )[, gene_name]
)


## generate the directory that will contain the figures
dir.create(snakemake@output[["methylation_figure_correlation_dir"]])


##########
## load the kegg metabolic pathways
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

## rename the kegg gene list with the kegg pathway names
names(KEGG_genes_list) <- KEGG_genes_file_names



###########
## create a data table that contain the KEGG information
###########

## extract all the gene associated with the KEGG
kegg_gene_vector <- unique(unlist(KEGG_genes_list))

## retreive the neurotransmission genes that we used previously
kegg_gene_vector <- intersect(kegg_gene_vector, neurotransmission_gene_vector)

## create the data table
kegg_information <- data.table(
    gene_name = kegg_gene_vector
)

## add the kegg pathway information into the data table kegg_information
for (i_kegg_list in seq(1, length(KEGG_genes_list))) {

    ## retreinve the kegg pathway name
    kegg_name <- names(KEGG_genes_list[i_kegg_list])

    ## extract the genes associated with the kegg pathway
    gene_of_interest <- unlist(KEGG_genes_list[i_kegg_list])

    ## update the kegg_information data table
    kegg_information[gene_name %in% gene_of_interest, kegg_pathway := kegg_name]

}


#############################################
## Plot of the average correlation between methylation and expression
#############################################


##########
## Formate the methylation data
##########


## melt the data to put all the probes into a same column
methylation_data_dt <- melt(
    methylation_data_dt,
    id.vars = c('sample_id', 'cluster'),
    measure.vars = colnames(methylation_data_dt)[!(colnames(methylation_data_dt) %in% c("sample_id", "cluster"))]
)

## do the average of methylation data per sample
methylation_data_average_dt <- methylation_data_dt[, .(methylation_value = mean(as.numeric(value), na.rm = TRUE)), by = c("sample_id", 'cluster')]

## center the methylation value in 0
#methylation_data_average_dt[, methylation_value := sapply(methylation_value, function(x) x - 0.5)]


##########
## Formate the expression data
##########


## transpose the expression data table
expression_data_dt <- transpose_datatable(
    expression_data_dt,
    column_name = 'Gene_Name',
    new_name = "sample_id"
)

## melt the data to put all the probes into a same column
expression_data_dt <- melt(
    expression_data_dt,
    id.vars = c('sample_id'),
    measure.vars = colnames(expression_data_dt)[!(colnames(expression_data_dt) %in% c("sample_id"))]
)

## do the average of expression data per sample
expression_data_average_dt <- expression_data_dt[, .(expression_value = mean(as.numeric(value), na.rm = TRUE)), by = "sample_id"]

## merge the data
merge_data_dt <- merge(
    methylation_data_average_dt,
    expression_data_average_dt,
    by = "sample_id",
    all.x = T
)


#############################################
## Correlation between methylation and expression per samples
#############################################

## extract the cluster name and put them into a vector
cluster_name_vector <- unique(unlist(merge_data_dt[, cluster]))

## initialize the data table that will contail the pvalue and the pearson coefficient
stats_corr_per_sample_dt <- data.table(
    cluster_name = character(),
    method_name = character(),
    coef = numeric(),
    pvalue = numeric()
)

## calculate the correlation within each clusteer
for (i_cluster in seq(1, length(cluster_name_vector))) {
    
    ## retrieve the cluster name
    cluster_name <- cluster_name_vector[i_cluster]

    ## extract the data associated with the cluster name
    data_subset <- copy(merge_data_dt[cluster == cluster_name,])

    ## calculate the pearson correlation and the pvalue
    corr_test <- cor.test(
        as.numeric(unlist(data_subset[,methylation_value])),
        as.numeric(unlist(data_subset[,expression_value])),
        use = "complete.obs",
        method = "pearson"
    )

    ## extract the pval 
    corr_pval <- corr_test$p.value

    ## extract the statistic
    corr_coef <- corr_test$estimate

    ## extract the method name
    corr_method <- corr_test$method

    ## create the table that contain the information of the cluster
    stats <- data.table(
        cluster_name = cluster_name,
        method_name = corr_method,
        coef = corr_coef,
        pvalue = corr_pval
    )

    ## merge the stat with the table that contain all the stats
    stats_corr_per_sample_dt <- rbind(
        stats_corr_per_sample_dt,
        stats
    )


}

## generate the average correlation plot
average_correlation_plot <- ggplot(
    data = merge_data_dt,
    aes(
        x = expression_value,
        y = methylation_value,
        fill = cluster
    )
) +
    geom_point(size = 4, pch = 21, color = 'black', stroke = 2,  alpha = 0.8) +
    #stat_ellipse(size = 3, aes(color = cluster)) +
    stat_smooth(
        aes(color = cluster),
        method = "lm",
        #col = "black",
        size = 3,
        formula = y ~ x,
        se = FALSE
    ) + 
    #geom_hline(yintercept = 0, color = "black") +
    #geom_vline(xintercept = 10, color = "black", size = 1) +
    scale_y_continuous(limits=c(0.35, 0.62)) +
    labs(
        x = "Average Gene Expression Value",
        y = "Average Gene Methylation Value"
    ) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    theme_bw(base_size = 24)


## save the plot
ggsave(
    filename = paste(
        snakemake@output[['methylation_figure_correlation_dir']],
        '/',
        "average_correlation_per_sample",
        '.svg',
        sep = ""
    ),
    plot = average_correlation_plot,
    device = "svg",
    width = 12,
    height = 10
)


#############################################
## Plot of the correlation between methylation and expression per genes
#############################################


## function for retreiving the gene name vector corresponding to the probe id vector 
get_genename_from_probeid <- function(
    probe_name_input
) {
    
    ## retreive the gene name associated with the probe id input
    gene_name_output <- methylation_information_dt[probe_name %in% probe_name_input]

    ## retorder the colum based on the input
    gene_name_output <- unlist(gene_name_output[match(probe_name_input, probe_name), gene_name])

    ## return the gene_name
    return(gene_name_output)
}

## put the gene name into the methylation data
methylation_data_dt[, variable := get_genename_from_probeid(variable)]

## do the average of methylation data per gene
methylation_data_average_gene_dt <- methylation_data_dt[, .(methylation_value = mean(as.numeric(value), na.rm = TRUE)), by = "variable"]

## do the average of expression data per gene
expression_data_average_gene_dt <- expression_data_dt[, .(expression_value = mean(as.numeric(value), na.rm = TRUE)), by = "variable"]

## merge the data
merge_data_gene_dt <- merge(
    methylation_data_average_gene_dt,
    expression_data_average_gene_dt,
    by = "variable",
    all.x = T
)

## merge the expression and methylation data with the kegg information
merge_data_gene_dt <- merge(
    merge_data_gene_dt,
    kegg_information,
    by.x = "variable",
    by.y = "gene_name",
    all.x = T
)

## calculate the pearson correlation and the pvalue
corr_test <- cor.test(
    as.numeric(unlist(merge_data_gene_dt[,methylation_value])),
    as.numeric(unlist(merge_data_gene_dt[,expression_value])),
    use = "complete.obs",
    method = "pearson"
)

## extract the pval 
corr_pval <- corr_test$p.value

## extract the statistic
corr_coef <- corr_test$estimate

## extract the method name
corr_method <- corr_test$method

## generate the title of the graphe with the corr test information
figure_name <- paste(
    corr_method,
    "\n",
    "r = ",
    corr_coef,
    "\n",
    "pvalue = ",
    corr_pval,
    sep = ""
)


## generate the average correlation plot
average_correlation_gene_plot <- ggplot(
    data = merge_data_gene_dt,
    aes(
        x = expression_value,
        y = methylation_value
        #fill = kegg_pathway
    )
) +
    geom_point(size = 4, pch = 21, color = 'black', fill = '#ebebeb', stroke = 2,  alpha = 0.8) +
    stat_smooth(size = 2, method = "lm", col = "red", formula = y ~ x, se = FALSE) + 
    labs(
        x = "Average Gene Expression Value",
        y = "Average Gene Methylation Value",
        title = "",
        subtitle = figure_name,
    ) +
    scale_color_manual(values = brewer.pal(9, "Set2")) +
    theme_bw(base_size = 24) +
    theme(legend.key.size = unit(1.5, 'cm'))


## save the plot
ggsave(
    filename = paste(
        snakemake@output[['methylation_figure_correlation_dir']],
        '/',
        "average_correlation_per_gene",
        '.svg',
        sep = ""
    ),
    plot = average_correlation_gene_plot,
    device = "svg",
    width = 13,
    height = 10
)


##########
## For each kegg metabolic pathways, create the correlation data 
##########


## initalize the list that will contain the kegg plots
kegg_plot_list <- list()

## Extract the probes associated with a specific metabolic pathway, from the kegg gene list
for (i_kegg in seq(1, length(KEGG_genes_list))) {

    ## retreive the name of the kegg pathway
    kegg_name <- names(KEGG_genes_list[i_kegg])

    ## extract the genes associated with a kegg pathway sublist
    kegg_gene_vector <- unlist(KEGG_genes_list[[i_kegg]])

    ## extract the methylation data associated with the kegg gene
    subset_merge_data <- merge_data_gene_dt[variable %in% kegg_gene_vector, ]

    ## calculate the pearson correlation and the pvalue
    corr_test <- cor.test(
        as.numeric(unlist(subset_merge_data[,methylation_value])),
        as.numeric(unlist(subset_merge_data[,expression_value])),
        use = "complete.obs",
        method = "pearson"
    )

    ## extract the pval 
    corr_pval <- corr_test$p.value

    ## extract the statistic
    corr_coef <- corr_test$estimate

    ## extract the method name
    corr_method <- corr_test$method

    ## generate the title of the graphe with the corr test information
    figure_name <- paste(
        corr_method,
        "\n",
        "r = ",
        corr_coef,
        "; ",
        "pvalue = ",
        corr_pval,
        sep = ""
    )

    ## generate the average correlation plot
    subset_merge_data <- ggplot(
        data = subset_merge_data,
        aes(
            x = expression_value,
            y = methylation_value
            #color = cluster
        )
    ) +
        geom_point(size = 5) +
        stat_smooth(method = "lm", col = "red", formula = y ~ x, se = FALSE) + 
        #geom_hline(yintercept = 0, color = "black", size = 1) +
        #scale_y_continuous(limits=c(-0.15, 0.15)) +
        labs(
            x = "Average Gene Expression Value",
            y = "Average Gene Methylation Value",
            title = kegg_name,
            subtitle = figure_name,
        ) +
        scale_color_manual(values = color_palette) +
        theme_classic(base_size = 24)

    ## generate the path of the file
    path_file <- paste(
        snakemake@output[['methylation_figure_correlation_dir']],
        '/',
        "average_correlation_per_gene_",
        kegg_name,
        '.svg',
        sep = ""
    )

    ## save the plot
    ggsave(
        filename = path_file,
        plot = subset_merge_data,
        device = "svg",
        width = 12,
        height = 10
    )
}


#############################################
## write the data
#############################################

## rename the column variable
merge_data_gene_dt <- merge_data_gene_dt[, setnames(.SD, "variable", "gene_name")]

## write the file that contain the average methylation and expression data for each genes
fwrite(
    merge_data_gene_dt,
    snakemake@output[["methylation_expression_data"]],
    sep = ","
)

## write the stats concerning the correlation per samples
path_file <- paste(
    snakemake@output[['methylation_figure_correlation_dir']],
    '/',
    "stats_average_correlation_per_sample",
    '.csv',
    sep = ""
)
fwrite(
    stats_corr_per_sample_dt,
    path_file,
    sep = ","
)
