
#################################################
## Generate the expression profile for each cluster
#################################################

#################################################
## Load the libraries
#################################################

library(data.table)
library(tidyverse)
library(RColorBrewer)
library(egg)


#################################################
## load the data
#################################################

## load the expression data
expression_data <- fread(
    snakemake@input[["normalized_expression_data_expressionFiltered"]],
    sep = ","
)

## load the genes related to GABA glutamate and calcium pathway
gene_vector <- unique(
    fread(
        snakemake@input[["genes_highExpressed"]],
        sep = ","
    )[, gene_name]
)

## load the cluster data
clusters <- fread(
    snakemake@input[["cluster_group"]],
    sep = ","
)
print(clusters)

##########
## add the Healthy sample into the cluster data
##########


## retreive Healthy samples
healthy_sample_vector <- c(
    "TCGA-06-0680",
    "TCGA-06-0675",
    "TCGA-06-0678",
    "TCGA-06-AABW",
    "TCGA-06-0681"
)

## creathe data associated with the Healthy samples
healthy_cluster_data <- data.table(
    "sample_id" = healthy_sample_vector,
    "cluster" = 'Healthy'
)

## rbind the cluster data
clusters <- rbind(
    clusters,
    healthy_cluster_data
)


## load the cluster names into a vector
cluster_names <- unique(unlist(clusters[, cluster]))

## create the directory that will contain the expression profile plots
dir.create(snakemake@output[["neurotransmission_expression_profile_dir"]])

##########
## Load the kegg information
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

##########
## Load the genes that are differentially expressed between clusters
##########

## list the file path from the directory that contain the DE genes
DE_results_path_vector <- list.files(
    snakemake@input[["DE_genes_dir"]],
    pattern = "filtered",
    full.names = TRUE
)

## read each results
DE_results_list <- lapply(
    DE_results_path_vector,
    function(x) fread(
        x,
        sep = ","
    )

)

## list the file name from the directory that contain the DE genes
DE_results_name_vector <- list.files(
    snakemake@input[["DE_genes_dir"]],
    pattern = "filtered",
    full.names = FALSE
)

## formate the file name
DE_results_name_vector <- str_replace_all(
    DE_results_name_vector,
    pattern = "filtered_TCGA_IDHall_",
    replacement = ""
)
DE_results_name_vector <- str_replace_all(
    DE_results_name_vector,
    pattern = ".csv",
    replacement = ""
)

## rename the list
names(DE_results_list) <- DE_results_name_vector


#################################################
## Functions
#################################################

#####
## Function for extract and formate the expression data from a vector of genes
####
get_expression_data <- function(
    gene_vector,
    expression_datatable
) {

    ## extract from the expression_data, the expression associated with the genes of interest
    expression_data_filtered <- expression_datatable[Gene_Name %in% gene_vector, ]

    ## extract the header as a data table
    dt_header <- (as.data.table(as.list(colnames(expression_data_filtered))))

    ## rename the header
    dt_header <- dt_header[
        ,
        setnames(
            .SD,
            colnames(dt_header),
            colnames(expression_data_filtered)
        )
    ]

    ## rbind the expression_data_filtered and do the transposition of the expression data
    expression_data_filtered <- as.data.table(t(
        rbindlist(
            list(
                dt_header,
                expression_data_filtered
            ),
            fill = F
        )
    ))

    ## rename the column by the gene names
    expression_data_filtered <- expression_data_filtered[
        ,
        setnames(
            .SD,
            colnames(expression_data_filtered),
            unlist(expression_data_filtered[V1 == "Gene_Name",])
        )
    ][
        !(Gene_Name == "Gene_Name"), # remove the row corresponding to the headers
    ][
        ,
        setnames( # rename column containing the sample id 
            .SD,
            "Gene_Name",
            "sample_id"
        )
    ]


    ## formate the expression_data_filtered datatable structure
    expression_data_melted <- melt.data.table(
        expression_data_filtered,
        id.vars = "sample_id",
        measure.vars = colnames(expression_data_filtered)[!(colnames(expression_data_filtered) %in% c("sample_id"))]
    )

    ## transform the values to numeric
    expression_data_melted[, value := sapply(value, as.numeric)]

    ## merge the data with the cluster informations
    expression_data_melted <- merge(
        expression_data_melted,
        clusters,
        by = "sample_id",
        all.y = T
    )

    ## trnasform into factor
    expression_data_melted <- copy(expression_data_melted[, sample_id := sapply(sample_id, as.factor)])

    ## add the average profile
    expression_data_melted <- rbind(
        expression_data_melted,
        expression_data_melted[, list(value = mean(value)), by = variable][, cluster := "mean_profile"][, sample_id := "mean_profile"]
    )

    ## divide the value for each gene by the average value
    expression_data_melted <- expression_data_melted[, list(value = value/(value[cluster == "mean_profile"]), cluster = cluster, sample_id = sample_id), by = c("variable")]

    ## do the average profile of each cluster group for each gene
    expression_data_melted_average <- expression_data_melted[, list(value = mean(value)), by = c("cluster", "variable")]

    ######
    ## order the genes based on standard deviation
    ######

    ## extract the gene set
    gene_vector <- unique(expression_data_melted_average[, variable])

 
    ## calculate the standard deviation for each genes
    standard_deviation_dt <- expression_data_melted_average[
        ,
        list(standard_dev = sd(value)),
        by = 'variable' 
    ][order(standard_dev, decreasing = TRUE),]

    ## extract the gene order
    gene_order <- unlist(standard_deviation_dt[, variable])

    



    
    ## for ordering based on cluster expression
    test <- dcast(
        expression_data_melted_average,
        variable ~ cluster,
        value.var = "value"
    )
    #setorderv(test, c("NT-1", "NT-2"),  order = c(-1, 1))
    setorderv(test, c("Healthy", "NT-1"),  order = c(-1, -1))

    gene_order <- unlist(test[, variable])








    ## update the order
    expression_data_melted_average$variable <- factor(
        expression_data_melted_average$variable,
        levels = gene_order
    )

    ## return the formated expression data 
    return(expression_data_melted_average)
}


##########
## function for generate the average expression plot
##########

## function for the generation of the average expression profile
expression_profile_plot_average <- function(
    expression_data_input,
    cluster_name_vector_input,
    color_table_input,
    title_name = "",
    gene_name = TRUE
) {


    if (gene_name == TRUE) {
        ## generate the expression profile plot
        p <- ggplot(
            data = expression_data_input[cluster %in% cluster_name_vector_input, ],
            aes(
                x = variable,
                y = value,
                color = cluster
            )
        ) +
            theme_classic(base_size = 90) +
            geom_point(alpha = 0.3) +
            geom_line(aes(group = cluster), alpha = 0.7, lwd = 5) +
            geom_point(data = expression_data_input[cluster == "mean_profile", ], aes(x = variable, y = value, color = cluster), color = "black", size = 5, alpha = 0) +
            geom_line(data = expression_data_input[cluster == "mean_profile", ], aes(x = variable, y = value, color = cluster, group = cluster), lwd = 5, color = "black") +
            theme_bw(base_size = 50) +
            labs(
                subtitle = title_name,
                x = "Genes",
                y = "Normalized Expression"
            ) +
            # scale_y_continuous(limits=c(2.5, 18.5)) + # set the scale of the y axe
            #scale_y_continuous(limits = c(0.50, 2)) + # set the scale of the y axe
            theme(
                # axis.text.x=element_blank(),
                legend.key.size = unit(1.5, "cm")
            ) +
            scale_color_manual(
                values = color_table_input[cluster %in% cluster_name_vector_input, color] # set the color
            )
    }

    if (gene_name == FALSE) {
        ## generate the expression profile plot
        p <- ggplot(
            data = expression_data_input[cluster %in% cluster_name_vector_input, ],
            aes(
                x = variable,
                y = value,
                color = cluster
            )
        ) +
            theme_classic(base_size = 90) +
            geom_point(alpha = 0.3) +
            geom_line(aes(group = cluster), alpha = 0.7, lwd = 5) +
            geom_point(data = expression_data_input[cluster == "mean_profile", ], aes(x = variable, y = value, color = cluster), color = "black", size = 5, alpha = 0) +
            geom_line(data = expression_data_input[cluster == "mean_profile", ], aes(x = variable, y = value, color = cluster, group = cluster), lwd = 5, color = "black") +
            theme_bw(base_size = 60) +
            labs(
                subtitle = title_name,
                x = "Neurotransmission-related Genes",
                y = "Average Normalized Expression"
            ) +
            # scale_y_continuous(limits=c(2.5, 18.5)) + # set the scale of the y axe
            #scale_y_continuous(limits = c(0.50, 2)) + # set the scale of the y axe
            theme(
                axis.text.x=element_blank(),
                #axis.line.x=element_line(),
                #axis.line.y=element_line(),
                legend.key.size = unit(2, "cm")
            ) +
            scale_color_manual(
                values = color_table_input[cluster %in% cluster_name_vector_input, color] # set the color
            )
    }
    

    ## generate the file_path
    file_path <- paste(
        snakemake@output[['expression_profile_plot_directory']],
        "/average_expression_profil_",
        paste(cluster_name_vector_input, collapse = "_"), # if the vector contain more than two values
        ".jpeg",
        sep = ""

    )
    
    ## save the plot
    ggsave(
        plot = p,
        filename = file_path,
        device = "jpeg",
        height = 20,
        width = 70,
        limitsize = F
    )

    print("Plot saved!")

    return(p)

}


###################################################
## for each metabolic pathway, generate the expression profil plots for all the genes
###################################################

## initialize the list that will contain the plots
expression_plot_list <- list()

## initialize the list that will contain the expression data
expression_data_list <- list()

## for each metabolic pathway gene set contained in the KEGG gene list
for (i_gene_set in seq(1, length(KEGG_genes_list))) {
    
    ## retreive the metabolic names
    kegg_name <- names(KEGG_genes_list[i_gene_set])
    
    ## retreive the genes related to the kegg pathway
    kegg_gene_set <- unlist(KEGG_genes_list[i_gene_set], use.names = FALSE)

    ## do the intersection with the significant genes that we identified
    #kegg_gene_set_significant <- intersect(
    #    gene_vector,
    #    kegg_gene_set
    #)

    ## for each significant kegg genes, extract the expression data related to these genes
    kegg_gene_set_expression <- get_expression_data(
        kegg_gene_set,
        expression_data
    )

    ## put the expression subset data into a list for counting 
    expression_data_list[[kegg_name]] <- kegg_gene_set_expression

    ## set colors that will be used for the generation of the plots
    color_palette <- brewer.pal(9, "Set1")

    ## name of the clusters that we want to do the expression profile plot
    cluster_name_vector <- unique(unlist(clusters[,cluster]))

    ## create a datatable containing the association of colors for each cluster name
    color_data <- data.table(
        cluster = cluster_name_vector,
        color = color_palette[1:(length(cluster_name_vector))]
    )

    ## set the order of the cluster, putting the Healthy condition at the end
    cluster_name_vector_vector <- unique(unlist(kegg_gene_set_expression$cluster))
    cluster_name_vector_vector <- sort(cluster_name_vector_vector[cluster_name_vector_vector != "Healthy"])
    cluster_name_vector_vector <- c(cluster_name_vector_vector, "Healthy")

    ## set the order of the cluster
    kegg_gene_set_expression$cluster <- factor(
        kegg_gene_set_expression$cluster,
        levels = cluster_name_vector_vector
    )

    ## create the average expression profile plot
    plot_all_expression_profile_average <- expression_profile_plot_average(
        expression_data_input = kegg_gene_set_expression,
        cluster_name_vector_input = cluster_names,
        color_table_input = color_data,
        title_name = kegg_name,
        gene_name = FALSE
    )

    ## create the file path that will be used for the saving
    file_path <- paste(
        snakemake@output[['neurotransmission_expression_profile_dir']],
        '/',
        "average_expression_profil_",
        kegg_name,
        '.svg',
        sep = ""
    )

    ## save the expression profile average into a file
    svg(
        filename = file_path,
        height = 20,
        width = 40
    )
    print(plot_all_expression_profile_average)
    dev.off()

    ## add the expression plot into the list
    expression_plot_list <- append(expression_plot_list, list(plot_all_expression_profile_average))

}

## generate the plot that contain the different kegg average expression plot
merged_average_expression_plot <- ggarrange(
    plots = expression_plot_list,
    ncol = 1
)

## save the expression profile average into a file
svg(
    filename = paste(
        snakemake@output[['neurotransmission_expression_profile_dir']],
        '/',
        "average_expression_profil_all",
        '.svg',
        sep = ""
    ),
    height = 40,
    width = 40
)
print(merged_average_expression_plot)
dev.off()


###################################################
## for each metabolic pathway, generate the expression profil plots for the gene
## that are differentially expressed
###################################################

##########
## Extract the genes that are differentially expressed
##########

## set the pvalue and fold change cutoff
pval_cutoff <- as.numeric(snakemake@params[["pval_cutoff"]])
logfoldchange_cutoff <- as.numeric(snakemake@params[["logfoldchange_cutoff"]])

## extract the probes filtering with the cutoff values
DE_results_dt <- Reduce(
    x = DE_results_list,
    function(list1, list2) {
        
        ## retreive the data table
        dt1 <- as.data.table(list1)
        dt2 <- as.data.table(list2)

        ## filter by the pvalue cutoff
        dt1 <- dt1[padj < pval_cutoff,]
        dt2 <- dt2[padj < pval_cutoff,]

        ## filter by the log fold change cutoff
        dt1 <- dt1[(log2FoldChange < -logfoldchange_cutoff) | (log2FoldChange > logfoldchange_cutoff),]
        dt2 <- dt2[(log2FoldChange < -logfoldchange_cutoff) | (log2FoldChange > logfoldchange_cutoff),]

        ## merge the dt1 and dt2 together
        dt3 <- rbind(dt1, dt2)
        
        ## return output
        return(dt3)
    }
)


## retreive the differential methylated sites
DE_genes_vector <- unique(unlist(DE_results_dt[, Gene_Name]))

##########
## for each metabolic pathway, generate the expression profil plots for the DE genes
##########

## initialize the list that will contain the plots
expression_plot_DE_list <- list()

## initialize the list that will contain the expression data
expression_data_DE_list <- list()

## for each metabolic pathway gene set contained in the KEGG gene list
for (i_gene_set in seq(1, length(KEGG_genes_list))) {
    
    ## retreive the metabolic names
    kegg_name <- names(KEGG_genes_list[i_gene_set])
    
    ## retreive the genes related to the kegg pathway
    kegg_gene_set <- unlist(KEGG_genes_list[i_gene_set], use.names = FALSE)

    ## do the intersection with the significant genes that we identified
    kegg_gene_set_significant <- intersect(
        DE_genes_vector,
        kegg_gene_set
    )

    ## for each significant kegg genes, extract the expression data related to these genes
    kegg_gene_set_expression <- get_expression_data(
        kegg_gene_set_significant,
        expression_data
    )

    ## put the expression subset data into a list for counting 
    expression_data_DE_list[[kegg_name]] <- kegg_gene_set_expression

    ## set colors that will be used for the generation of the plots
    color_palette <- brewer.pal(9, "Set1")

    ## name of the clusters that we want to do the expression profile plot
    cluster_name_vector <- unique(unlist(clusters[,cluster]))

    ## create a datatable containing the association of colors for each cluster name
    color_data <- data.table(
        cluster = cluster_name_vector,
        color = color_palette[1:(length(cluster_name_vector))]
    )

    ## set the order of the cluster, putting the Healthy condition at the end
    cluster_name_vector_vector <- unique(unlist(kegg_gene_set_expression$cluster))
    cluster_name_vector_vector <- sort(cluster_name_vector_vector[cluster_name_vector_vector != "Healthy"])
    cluster_name_vector_vector <- c(cluster_name_vector_vector, "Healthy")

    ## set the order of the cluster
    kegg_gene_set_expression$cluster <- factor(
        kegg_gene_set_expression$cluster,
        levels = cluster_name_vector_vector
    )


    ## generate the title of the figure
    figure_title <- kegg_name
    figure_title <- str_replace(figure_title, pattern = "_", replacement = " ")

    ## create the average expression profile plot
    plot_all_expression_profile_average <- expression_profile_plot_average(
        expression_data_input = kegg_gene_set_expression,
        cluster_name_vector_input = cluster_names,
        color_table_input = color_data,
        title_name = figure_title,
        gene_name = FALSE
    )

    ## create the file path that will be used for the saving
    file_path <- paste(
        snakemake@output[['neurotransmission_expression_profile_dir']],
        '/',
        "DE_average_expression_profil_",
        kegg_name,
        '.svg',
        sep = ""
    )

    ## save the expression profile average into a file
    svg(
        filename = file_path,
        height = 20,
        width = 40
    )
    print(plot_all_expression_profile_average)
    dev.off()

    ## add the expression plot into the list
    expression_plot_DE_list <- append(expression_plot_DE_list, list(plot_all_expression_profile_average))

}

### generate the plot that contain the different kegg average expression plot
#merged_average_expression_plot <- ggarrange(
#    plots = expression_plot_DE_list,
#    ncol = 1
#)

### save the expression profile average into a file
#svg(
#    filename = paste(
#        snakemake@output[['neurotransmission_expression_profile_dir']],
#        '/',
#        "DE_average_expression_profil_all",
#        '.svg',
#        sep = ""
#    ),
#    height = 40,
#    width = 40
#)
#print(merged_average_expression_plot)
#dev.off()


#################################################
## Generate the count for each metabolic pathway
#################################################


## function for counting the number of time the cluster is associated with the higher or the lower gene expression
counting_genes <- function(
    average_expression_data, # output generated wi the get_expression_data function
    pathway_name # name of the pathway that is associated with the expression data
) {

    ## remove the lines related to the healthy samples and average
    #average_expression_data <- average_expression_data[!(cluster %in% c("Healthy", "mean_profile")),]



    #####
    ## Do the counting for healthy samples
    #####
    
    average_expression_data <- average_expression_data[!(cluster %in% c("mean_profile")),]







    ## initialize the data table that will contain the cluster names
    cluster_name_dt <- data.table(
        cluster_name = unique(unlist(unlist(average_expression_data[,cluster])))
    )

    ## initialize the list that will contain the counting list output
    counting_output_list <- list()

    ## retreive the unique vector that contain the gene name from the average expression data
    gene_name_vector <- unique(unlist(average_expression_data[, variable]))

    ## number of gene in total
    gene_total <- length(gene_name_vector)

    ## initalize the table that contain the counting
    counting_output_dt <- data.table(
        "gene_name" = character(),
        "min_expression_cluster" = character(), # cluster name associated with the lower expression value
        "max_expression_cluster" = character() # cluster name associated with the higher expression value
    )

    ## extract the data associated with each gene from the gene_name_vector
    for (i_gene in seq(1, length(gene_name_vector))) {

        ## retreive the name of the gene
        gene_name <- gene_name_vector[i_gene]

        ## extract subset expression data
        subset_expression_data <- average_expression_data[variable == gene_name, ]

        ## order the subset data with the expression value
        subset_expression_data <- subset_expression_data[order(value),]

        ## put the cluster names with higher and lower expression value
        counting_dt <- data.table(
            "gene_name" = gene_name,
            "min_expression_cluster" = unlist(subset_expression_data[1, cluster]), 
            "max_expression_cluster" = unlist(subset_expression_data[nrow(subset_expression_data), cluster])
        )

        ## rbind the counting_dt with the counting_output_dt
        counting_output_dt <- rbind(
            counting_output_dt,
            counting_dt
        )
    }

    ## do the counting
    counting_min <- as.data.table(table(counting_output_dt[,min_expression_cluster]))[, setnames(.SD, "V1", "cluster_name")]
    counting_max <- as.data.table(table(counting_output_dt[,max_expression_cluster]))[, setnames(.SD, "V1", "cluster_name")]

    ## merge with the cluster data with the counting data
    counting_min <- merge(
        cluster_name_dt,
        counting_min,
        by = "cluster_name",
        all.x = T
    )
    counting_max <- merge(
        cluster_name_dt,
        counting_max,
        by = "cluster_name",
        all.x = T
    )
    
    ## put a zero value when there is NA
    counting_min[is.na(counting_min)] <- 0
    counting_max[is.na(counting_max)] <- 0

    ## do the ratio
    counting_min[, ratio := paste(N, "/", gene_total, sep = "")]
    counting_max[, ratio := paste(N, "/", gene_total, sep = "")]

    ## calculate the percentage
    counting_min[, percentage := round(N/gene_total*100, digits = 2)]
    counting_max[, percentage := round(N/gene_total*100, digits = 2)]

    ## put the counting max and min into the output list
    counting_output_list[["min"]] <- (counting_min)
    counting_output_list[["max"]] <- (counting_max)

    ## put min and max in the data tables
    counting_min[, type := "underexpression"]
    counting_max[, type := "overexpression"]

    ## rbind the min and max data
    counting_output_dt <- rbind(
        counting_min,
        counting_max
    )

    ## put the name the pathway
    counting_output_dt[, pathway := pathway_name]

    ## return the list output
    return(counting_output_dt)

}

##########
## count for each pathway data expression
##########

## initialize the data table that will contain all the data
counting_pathway_dt <- data.table(
    "cluster_name" = character(),
    "N" = numeric(),
    "ratio" = character(),
    "percentage" = numeric(),
    "type" = character(),
    "pathway" = character()
)

for (i_expression_data in seq(1, length(expression_data_list))) {

    ## retreive the pathway name
    pathway_name <- names(expression_data_list[i_expression_data])

    ## retreive the expression data related to the pathway_name
    pathway_expression_data <- as.data.table(expression_data_list[i_expression_data])

    ## rename the expression data column
    colnames(pathway_expression_data) <- c(
        "cluster",
        "variable",
        "value"
    )

    ## do the counting of the expression data related to the pathway
    pathway_counting <- counting_genes(pathway_expression_data, pathway_name)

    ## merge with the data together
    counting_pathway_dt <- rbind(
        counting_pathway_dt,
        pathway_counting
    )

}


##########
## count for each pathway data expression for the DE genes
##########


## initialize the data table that will contain all the data
counting_pathway_DE_dt <- data.table(
    "cluster_name" = character(),
    "N" = numeric(),
    "ratio" = character(),
    "percentage" = numeric(),
    "type" = character(),
    "pathway" = character()
)



for (i_expression_data in seq(1, length(expression_data_DE_list))) {

    ## retreive the pathway name
    pathway_name <- names(expression_data_DE_list[i_expression_data])

    ## retreive the expression data related to the pathway_name
    pathway_expression_data <- as.data.table(expression_data_DE_list[i_expression_data])
    
    ## rename the expression data column
    colnames(pathway_expression_data) <- c(
        "cluster",
        "variable",
        "value"
    )

    #print("------------")
    #print(pathway_name)
    #print(pathway_expression_data)
    #print(table(pathway_expression_data[,variable]))

    ## do the counting of the expression data related to the pathway
    pathway_counting <- counting_genes(pathway_expression_data, pathway_name)

    ## merge with the data together
    counting_pathway_DE_dt <- rbind(
        counting_pathway_DE_dt,
        pathway_counting
    )

}


##########
## Write the other data
##########


## write the counting data table into file
file_path <- paste(
    snakemake@output[['neurotransmission_expression_profile_dir']],
    '/',
    "counting_data",
    '.csv',
    sep = ""
)
fwrite(
    counting_pathway_dt,
    file_path,
    sep = ","
)

file_path <- paste(
    snakemake@output[['neurotransmission_expression_profile_dir']],
    '/',
    "DE_counting_data",
    '.csv',
    sep = ""
)
fwrite(
    counting_pathway_DE_dt,
    file_path,
    sep = ","
)

#Sys.sleep(1000000)


