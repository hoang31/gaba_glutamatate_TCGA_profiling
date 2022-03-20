
#################################################
## Generate some figures from the methylation data
#################################################


#################################################
## Load the libraries
#################################################


library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(egg)


## set the color palette
color_palette <- brewer.pal(9, "Set1")

#################################################
## Load the data 
#################################################


## generate the directory that will contain the figures
dir.create(snakemake@output[["methylation_figure_dir"]])

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

## read the data table that contain the methylation information
methylation_information <- fread(
    snakemake@input[["methylation_information"]],
    sep = ","
)

## modification of the methylation information
methylation_information[grep(gene_name, pattern = 'SLC25A6'), gene_name := substr(gene_name, start = 1, stop = 7)]


## extract all the file path from the methylation data directory
methylation_data_path_vector <- list.files(
    snakemake@input[["DMgenes_analysis_dir"]],
    full.names = T
)

## retreive the methylation data
methylation_data_dt <- fread(
    grep(methylation_data_path_vector, pattern = "methylation_data_filtered.csv", value = T),
    sep = ","
)

## rename the cluster name for healthy samples
healthy_sample_vector <- c(
    "TCGA-06-AABW",
    "TCGA-74-6573"
)

methylation_data_dt[sample_id %in% healthy_sample_vector, cluster := "Healthy"]


##########
## Retreive the differential methylated sites
##########


## retreive the differential methylation data
diff_methylation_data_list <- lapply(
    grep(methylation_data_path_vector, pattern = "_vs_", value = T),
    function(x) fread(
        x,
        sep = ","
    )

)

## extract all the file name from the methylation data directory
methylation_data_name_vector <- list.files(
    snakemake@input[["DMgenes_analysis_dir"]],
    full.names = FALSE,
    pattern = "_vs_"
)

## rename the sublist
names(diff_methylation_data_list) <- methylation_data_name_vector

## remove the .csv at the end of the name
names(diff_methylation_data_list) <- str_replace(
    names(diff_methylation_data_list),
    pattern = ".csv",
    replacement = ""
)

## load the genes related to GABA glutamate and calcium pathway
neurotransmission_gene_vector <- unique(
    fread(
        snakemake@input[["genes_highExpressed"]],
        sep = ","
    )[, gene_name]
)

#################################################
## ANALYSIS
#################################################


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

###########
## Extract the probes that are significant
###########


## set the pvalue and fold change cutoff
pval_cutoff <- as.numeric(snakemake@params[["pval_cutoff"]])
logfoldchange_cutoff <- as.numeric(snakemake@params[["logfoldchange_cutoff"]])

## extract the probes filtering with the cutoff values
diff_methylated_dt <- Reduce(
    x = diff_methylation_data_list,
    function(list1, list2) {
        
        ## retreive the data table
        dt1 <- as.data.table(list1)
        dt2 <- as.data.table(list2)

        ## filter by the pvalue cutoff
        dt1 <- dt1[adj_pvalue < pval_cutoff,]
        dt2 <- dt2[adj_pvalue < pval_cutoff,]

        ## filter by the log fold change cutoff
        dt1 <- dt1[(logfoldchange < -logfoldchange_cutoff) | (logfoldchange > logfoldchange_cutoff),]
        dt2 <- dt2[(logfoldchange < -logfoldchange_cutoff) | (logfoldchange > logfoldchange_cutoff),]

        ## merge the dt1 and dt2 together
        dt3 <- rbind(dt1, dt2)[order(logfoldchange),]
        
        ## return output
        return(dt3)
    }
)

## retreive the differential methylated sites
diff_methylated_probe_vector <- unique(unlist(diff_methylated_dt[, probe_id]))

## extract the methylation data associated with the DM probes
methylation_data_dt_differential <- copy(methylation_data_dt[, c("sample_id", "cluster", diff_methylated_probe_vector), with = F])

# retreive the probes contained from the methylation data
probe_name_of_interest_vector <- colnames(methylation_data_dt_differential)[!(colnames(methylation_data_dt_differential) %in% c("sample_id", "cluster"))]


###########
## calculate the average methylation
###########


# do the average in the same column all the values
methylation_data_melted <- melt(
    methylation_data_dt_differential,
    id.vars = 'cluster',
    measure.vars = colnames(methylation_data_dt_differential)[!(colnames(methylation_data_dt_differential) %in% c("sample_id", "cluster"))]
)


## function for retreiving the gene name vector corresponding to the probe id vector 
get_genename_from_probeid <- function(
    probe_name_input
) {
    
    ## retreive the gene name associated with the probe id input
    gene_name_output <- methylation_information[probe_name %in% probe_name_input]

    ## retorder the colum based on the input
    gene_name_output <- unlist(gene_name_output[match(probe_name_input, probe_name), gene_name])

    ## return the gene_name
    return(gene_name_output)
}

## put the gene name into the methylation data
methylation_data_melted[, gene_name := get_genename_from_probeid(variable)]

## function for retreiving the kegg pathway from the gene name vector
get_kegg_from_genename <- function(
    gene_name_vector
) {
    
    ## retreive the gene name associated with the probe id input
    kegg_vector_output <- kegg_information[gene_name %in% gene_name_vector]

    ## reorder the colum based on the input
    kegg_vector_output <- unlist(kegg_vector_output[match(gene_name_vector, gene_name), kegg_pathway])

    ## return the gene_name
    return(kegg_vector_output)
}

## put the kegg pathway information into the methylation data
methylation_data_melted[, kegg_pathway := get_kegg_from_genename(gene_name)]

## do the average by cluster
methylation_data_average <- methylation_data_melted[, .(average = mean(as.numeric(value), na.rm = TRUE)), by = c("cluster", "kegg_pathway")]

## rename pathways
methylation_data_average[, kegg_pathway := str_replace(kegg_pathway, pattern = "[_]", replacement = "\n")]
methylation_data_melted[, kegg_pathway := str_replace(kegg_pathway, pattern = "[_]", replacement = "\n")]

## center the methylation in 0
#methylation_data_average[, average := sapply(average, function(x) x - 0.5)]


#################################################
## do the statistics for each gene
#################################################


## calculate the methylation value per gene and cluster
methylation_data_per_gene_per_cluster <- methylation_data_melted[, .("value" = mean(value, na.rm = T)), by = c("gene_name", "cluster")]

## retreive the cluster names
cluster_name_vector <- unique(unlist(methylation_data_per_gene_per_cluster[, cluster]))

## do the combination of the cluster name vector
cluster_combination <- combn(
    cluster_name_vector,
    m = 2
)

## initialize the data table that will contain the pvalue
stats_dt <- data.table(
    condition1 = character(),
    condition2 = character(),
    pvalue = numeric()
)

## for all combination, do the statistics
for (i_cluster_combination in seq(1, ncol(cluster_combination))) {

    ## retreive the combination
    condition1 <- cluster_combination[1, i_cluster_combination]
    condition2 <- cluster_combination[2, i_cluster_combination]

    ## retreive the data associated with the conditions
    data_condition1 <- unlist(methylation_data_per_gene_per_cluster[cluster == condition1, value])
    data_condition2 <- unlist(methylation_data_per_gene_per_cluster[cluster == condition2, value])
    
    ## do the stats and extract the pvalue
    pval <- wilcox.test(
        data_condition1,
        data_condition2
    )$p.value

    ## generate a table that contain the information of the stat
    stats_subset_dt <- data.table(
        condition1 = condition1,
        condition2 = condition2,
        pvalue = pval
    )

    ## merge the table with the stats_dt
    stats_dt <- rbind(
        stats_dt,
        stats_subset_dt
    )
}

## calculate the pvalue corrected by bonferonni
stats_dt[, adj_pvalue := p.adjust(pvalue, method = "bonferroni", n = nrow(stats_dt))]

## generate the file path of the file that will contain the pvalues
stats_file_path <- paste(
    snakemake@output[["methylation_figure_dir"]],
    "/",
    "boxplot_pvalue",
    ".csv",
    sep = ""
)

## write the data table into a file
fwrite(
    stats_dt,
    stats_file_path,
    sep = ","
)

#################################################
## calculate the average methylation value for each genes by metabolic pathways
#################################################


## melth the methylation data
methylation_data_melted2 <- melt(
    methylation_data_dt_differential,
    id.vars = 'cluster',
    measure.vars = colnames(methylation_data_dt_differential)[!(colnames(methylation_data_dt_differential) %in% c("sample_id", "cluster"))]
)

## do the average of methylation value by probes and by clusters
methylation_data_average_by_probe <- methylation_data_melted2[, .(average = mean(as.numeric(value), na.rm = TRUE)), by = c("variable", "cluster")][order(variable),]

## center the methylation in 0
#methylation_data_average_by_probe[, average := sapply(average, function(x) x - 0.5)]

## extract the gene name associated with the probes of intrest
methylation_data_average_by_probe[, gene_name := get_genename_from_probeid(variable)]

## initalize the list that will contain the kegg plots
kegg_plot_list <- list()

## Extract the probes associated with a specific metabolic pathway, from the kegg gene list
for (i_kegg in seq(1, length(KEGG_genes_list))) {

    ## retreive the name of the kegg pathway
    kegg_name <- names(KEGG_genes_list[i_kegg])

    ## extract the genes associated with a kegg pathway sublist
    kegg_gene_vector <- unlist(KEGG_genes_list[[i_kegg]])

    ## retreive the gene of interest that we used previously
    kegg_gene_vector <- intersect(neurotransmission_gene_vector, kegg_gene_vector)

    ## extract the methylation data associated with the kegg gene
    kegg_probes_methylation_data <- methylation_data_average_by_probe[gene_name %in% kegg_gene_vector, ]

    ## do the average by gene and cluster
    meth_average_data <- kegg_probes_methylation_data[, .(average = mean(as.numeric(average), na.rm = TRUE)), by = c("gene_name", "cluster")]

    ## set the order of the genes for the vizualization
    standard_deviation_dt <- meth_average_data[
            ,
            list(standard_dev = sd(average)),
            by = 'gene_name' 
        ][order(standard_dev, decreasing = TRUE),]

    

    ## extract the gene order
    gene_order <- unlist(standard_deviation_dt[, gene_name])

    ### dcast the data for ordering by the NT-1 and NT-2 cluster
    #test <- dcast(
    #    meth_average_data,
    #    gene_name ~ cluster,
    #    value.var = "average"
    #)
    #setorderv(test, c("NT-2", "NT-1"))
    #gene_order <- unlist(test[, gene_name])

    ## order by the NT-3 cluster
    gene_order <- unlist(meth_average_data[cluster == "NT-3",][order(average, decreasing = T),gene_name])

    ## update the order
    meth_average_data$gene_name <- factor(
        meth_average_data$gene_name,
        levels = gene_order
    )

    ## set the order of the cluster, putting the Healthy condition at the end
    cluster_name_vector_ordered <- unique(unlist(meth_average_data$cluster))
    cluster_name_vector_ordered <- sort(cluster_name_vector_ordered[cluster_name_vector_ordered != "Healthy"])
    cluster_name_vector_ordered <- c(cluster_name_vector_ordered, "Healthy")

    ## set the order of the cluster
    meth_average_data$cluster <- factor(
        meth_average_data$cluster,
        levels = cluster_name_vector_ordered
    )

    ## generate the title of the plot
    plot_title <- str_replace(kegg_name, pattern = "_", replacement = " ")

    ## generate the figure associated with the methylation value by probes
    methylation_profiles <- ggplot(
        data = meth_average_data,
        aes(
            x = gene_name,
            y = average,
            color = cluster
        )
    ) +
        #geom_bar(stat="identity", position = position_dodge()) +
        geom_point(alpha = 0.3) +
        geom_line(aes(group = cluster), alpha = 0.7, lwd = 5) +
        geom_hline(yintercept = 0, color = "black") +
        scale_y_continuous(limits=c(0, 1)) +
        scale_color_manual(values = color_palette) +
        labs(
            title = plot_title,
            x = "Neurotransmission-related Genes",
            y = "Average Methylation Levels\n(Beta Values)"
        ) +
        theme_bw(base_size = 60) +
        theme(
            axis.text.x = element_blank(),
            legend.key.size = unit(2, 'cm')

        )
    
    ## generate the path of the file
    file_path <- paste(
        snakemake@output[["methylation_figure_dir"]],
        "/",
        "methylation_profiles_",
        kegg_name,
        ".svg",
        sep = ""
    )

    ## save the file into a figure
    ggsave(
        methylation_profiles,
        filename = file_path,
        device = "svg",
        height = 8,
        width = 40,
        bg = "white"
    )

    ## add the expression plot into the list
    kegg_plot_list <- append(
        kegg_plot_list,
        list(methylation_profiles)
    )

}



## generate the plot that contain the different kegg methylation plot
methylation_profiles_all_kegg_plot <- ggarrange(
    plots = kegg_plot_list,
    ncol = 1
)

## save the methylation profile average into a file
svg(
    filename = paste(
        snakemake@output[['methylation_figure_dir']],
        '/',
        "methylation_profiles_all_kegg",
        '.svg',
        sep = ""
    ),
    height = 40,
    width = 40
)
print(methylation_profiles_all_kegg_plot)
print(methylation_profiles_all_kegg_plot)

dev.off()


#################################################
## GENERATER THE FIGURES
#################################################

##########
## barchart of methylation data
##########

## generate the plot of the average methylation data
methylation_barchart_average <- ggplot(
    data = methylation_data_average,
    aes(
        x = kegg_pathway,
        y = as.numeric(average),
        fill = cluster
    )
) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_y_continuous(limits=c(0, 1)) +
    scale_fill_manual(values = color_palette) +
    labs(
        x = "Kegg Metabolic Pathways",
        y = "Average Methylation Levels\n(Beta Values)"
    ) +
    theme_bw(base_size = 24)

ggsave(
    methylation_barchart_average,
    filename = paste(
        snakemake@output[['methylation_figure_dir']],
        '/',
        "methylation_average_barchart",
        '.svg',
        sep = ""
    ),
    device = "svg",
    height = 8,
    width = 18
)

##########
## boxplot of methylation data
##########


## set the order of the cluster, putting the Healthy condition at the end
cluster_name_vector_ordered <- unique(unlist(methylation_data_melted$cluster))
cluster_name_vector_ordered <- sort(cluster_name_vector_ordered[cluster_name_vector_ordered != "Healthy"])
cluster_name_vector_ordered <- c(cluster_name_vector_ordered, "Healthy")

## set the order of the cluster
methylation_data_melted$cluster <- factor(
    methylation_data_melted$cluster,
    levels = cluster_name_vector_ordered
)


### generate the boxplots of all the methylation data
methylation_boxplot <- ggplot(
    data = methylation_data_melted,
    aes(
        x = kegg_pathway,
        y = as.numeric(value),
        fill = cluster
    )
) +
    geom_boxplot(size = 2, color = 'black', notch = T) +
    #geom_jitter(color="black", size=0.4, alpha=0.5) +
    scale_fill_manual(values = color_palette) +
    labs(
        x = "Kegg Metabolic Pathways",
        y = "Methylation Levels\n(Beta Values)"
    ) +
    theme_bw(
        base_size = 24
    ) + 
    theme(
        legend.key.size = unit(1.5, 'cm')
    )


ggsave(
    methylation_boxplot,
    filename = paste(
        snakemake@output[['methylation_figure_dir']],
        '/',
        "methylation_average_boxplot",
        '.svg',
        sep = ""
    ),
    device = "svg",
    height = 8,
    width = 24
)


##########
## Generate a heatmap based on the methylation data
##########


## retreive the columns that are associated with only na values
na_position <- apply(
    methylation_data_dt_differential,
    2,
    function(x) sum(is.na(x))
) == nrow(methylation_data_dt_differential)

## filter the differential methylation data
methylation_data_filtered_dt <- methylation_data_dt_differential[, !na_position, with = F]

## transform the methylation data to matrix
methylation_matrix <- as.matrix(methylation_data_filtered_dt[, !c("cluster", "sample_id"), with = F])

## extract the information cluster associated with the samples
cluster_data <- methylation_data_dt[, c("cluster", "sample_id"), with = F]

## function for creating the annotation of the heatmap
create_annotation <- function(
    variable,
    clinicalData,
    horizontal_legend = F
) {
    ## create the colors
    all_colors <- c(
        brewer.pal(9, "Set1"),
        brewer.pal(8, "Dark2"),
        brewer.pal(8, "Set2"),
        brewer.pal(12, "Paired")
    )
    # print(all_colors)
    # print(all_colors)


    ## extract data table with the specific variables
    clinicalData <- clinicalData[, variable, with = F]

    ## initialise the color list
    colors_for_variables <- list()

    ## for each data associated with each variables, do a counting of occurence and extract colors
    for (i in (1:length(variable)))
    {

        ## extract the number of catagories in the variable i-th
        categories <- names(table(clinicalData[, variable[i], with = F]))

        ## take the i-th first colors from the "all_colors" variable
        colors_for_categories <- all_colors[1:length(categories)]

        ## remove the colors used from the initial color vector
        all_colors <- all_colors[!(all_colors %in% colors_for_categories)]
 

        ## name each color by the category
        names(colors_for_categories) <- categories

        ## append into the list the colors categories
        colors_for_variables[[i]] <- colors_for_categories
    }

    ## give variable names for each vector in the list
    names(colors_for_variables) <- variable

    ## create vertical annotation
    if (horizontal_legend == FALSE) {
        annotation <- HeatmapAnnotation(
            df = clinicalData,
            col = colors_for_variables
        )
    }

    ## create horizontal annotation
    if (horizontal_legend == TRUE) {
    ## create the annotations
        annotation <- HeatmapAnnotation(
            df = clinicalData,
            col = colors_for_variables,
            annotation_legend_param = list(
                nrow = 1
            )
        )
    }

    ## return the anntoation
    return(annotation)
}

## function for creating the heatmap based on the methylation data
create_heatmap_methylation <- function(
    methylation_matrix,
    clinical,
    variables,
    path_file,
    clustering_distance_rows = "minkowski",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "spearman",
    clustering_method_columns = "ward.D2"
) {

    ## generate the heatmap file
    generation_heatmap <- FALSE

    ## do the transpose of the matrix
    methylation_matrix <- t(methylation_matrix)
    
    if (!(missing(path_file))) {
        generation_heatmap <- TRUE
    }

    ## extract and make the annotations matrix
    annotation <- create_annotation(variables, clinical)

    ## generate the heatmap
    ht <- Heatmap(
        methylation_matrix,
        border = TRUE,
        #column_title = column_title,
        #row_title = row_title,
        #column_title_gp = column_title_gp,
        #row_title_gp = row_title_gp,
        show_column_names = FALSE,
        show_row_names = FALSE,
        
        ## legends for the heatmap
        top_annotation = annotation,
        heatmap_legend_param = list(
            title = "Methylation Beta Values",
            title_position = "leftcenter-rot",
            at = c(0, 0.5, 1),
            grid_width = unit(0.75, "cm"),
            legend_height = unit(6, "cm"),
            border = "black"
        ),

        ## number of clusters
        column_split = 4,
        column_gap = unit(2, "mm"),

        ## methods for clustering
        clustering_distance_rows = clustering_distance_rows,
        clustering_method_rows = clustering_method_rows,
        clustering_distance_columns = clustering_distance_columns,
        clustering_method_columns = clustering_method_columns,
    )
    
    ## write the figure into a file
    if (generation_heatmap) {
        svg(path_file, width = 12, height = 10) 
        draw(
            ht, 
            row_title = "Neurotransmission-related Genes",
            column_title = "Samples",
            column_title_gp = gpar(fontsize = 28, fontface = "bold"),
            row_title_gp = gpar(fontsize = 28, fontface = "bold")
        )
        dev.off()
    }

    ## return output
    return(ht)
}

## call the function for 
create_heatmap_methylation(
    methylation_matrix,
    cluster_data,
    c("cluster"),
    path_file = paste(
        snakemake@output[['methylation_figure_dir']],
        '/',
        "methylation_heatmap",
        '.svg',
        sep = ""
    ),
)




