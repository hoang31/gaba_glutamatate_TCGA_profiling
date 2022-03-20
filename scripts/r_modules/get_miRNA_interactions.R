

############################################
## Get the miRNA interactions with the neurotransmission related genes of interest
############################################


############################################
## load the libraries 
############################################


library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(egg)
#library(ggpubr)



## set the color palette
color_palette <- brewer.pal(9, "Set1")


############################################
## load the data 
############################################

## load the utils script
source(snakemake@params[["utils"]])

## load the data associated with the gene of interest
gene_of_interest_data_dt <- fread(
    snakemake@input[["methylation_expression_data"]],
    sep = ","
)

## load the data associated with the edge information
edge_information_dt <- fread(
    snakemake@input[["edge_information"]],
    sep = ","
)

## load the expression data related to the miRNA
miRNA_expression_data_dt <- fread(
    snakemake@input[["miRNA_expression_data"]],
    sep = ","
)


### load the expression data related to the miRNA of Healthy samples
#miRNA_expression_data_healthy_dt <- fread(
#    snakemake@input[["miRNA_expression_data_healthy"]],
#    sep = ","
#)

## load the cluster data
cluster_data_dt <- fread(
    snakemake@input[["cluster_group"]],
    sep = ","

)

## create the directory that will contain the figures
dir.create(snakemake@output[["miRNA_figure_dir"]])



##########
## load the deseq2 results
##########

## extract all the file path from the de results directory
mir_deseq_results_path_vector <- list.files(
    snakemake@input[["results_deseq2_miRNA"]],
    full.names = T
)

## retreive the differential methylation data
mir_deseq_results_dt <- lapply(
    mir_deseq_results_path_vector,
    function(x) fread(
        x,
        sep = ","
    )
)

## extract all the file name from the methylation data directory
mir_deseq_results_name_vector <- list.files(
    snakemake@input[["results_deseq2_miRNA"]],
    full.names = F
)

## rename the sublist
names(mir_deseq_results_dt) <- mir_deseq_results_name_vector

## remove the .csv at the end of the name
names(mir_deseq_results_dt) <- str_replace(
    names(mir_deseq_results_dt),
    pattern = ".csv",
    replacement = ""
)


############################################
## ANALYSIS
############################################


##########
## merge cancer and Healthy information data together
##########


### add the cluster information for the Healthy samples
#healthy_cluster_data <- data.table(
#    "sample_id" = colnames(miRNA_expression_data_healthy_dt)[(colnames(miRNA_expression_data_healthy_dt) != "gene_name")],
#    "cluster" = 'Healthy'
#)

### rbind the cluster data
#cluster_data_dt <- rbind(
#    cluster_data_dt,
#    healthy_cluster_data
#)

### merge the miRNA expression data
#miRNA_expression_data_dt <- merge(
#    miRNA_expression_data_dt,
#    miRNA_expression_data_healthy_dt,
#    by = "gene_name"
#)


############################################
## Get the miRNA of interest
############################################


## extract the genes of interest
gene_of_interest_vector <- unlist(gene_of_interest_data_dt[, gene_name])

## extract the edge informations related to the gene of interest
edge_information_of_interest_dt <- rbind(
    edge_information_dt[node1 %in% gene_of_interest_vector,],
    edge_information_dt[node2 %in% gene_of_interest_vector,]
)

## extract the edge informations related to miRNA
edge_information_of_interest_mir_dt <- rbind(
    edge_information_of_interest_dt[grep(node1, pattern = "miR|mir"),],
    edge_information_of_interest_dt[grep(node2, pattern = "miR|mir"),]
)

## count the number of interaction with a miRNA
mir_counting_dt <- edge_information_of_interest_mir_dt[, .(miRNA_interaction_count = .N), by = "node2"][,setnames(.SD, "node2", "gene_name")]

## merge the counting with the gene of interest information
gene_of_interest_data_dt <- merge(
    gene_of_interest_data_dt,
    mir_counting_dt,
    by = 'gene_name',
    all.x = T,
)[order(miRNA_interaction_count),]

## extract the miRNA of itnerest
mir_of_interest <- unique(
    c(
        unlist(edge_information_of_interest_mir_dt[grep(node1, pattern = "miR|mir"),node1]),
        unlist(edge_information_of_interest_mir_dt[grep(node2, pattern = "miR|mir"),node2])
    )
)

## formate the name of the miRNA gene name
#mir_of_interest <- str_replace_all(mir_of_interest, pattern = "-3p", "")
#mir_of_interest <- str_replace_all(mir_of_interest, pattern = "-5p", "")
mir_of_interest <- tolower(mir_of_interest)











#print(sort(mir_of_interest))
#break
## retreive the miRNA that does not have the -1 or -2 at the end of their
#mir_of_interest <- c(
#    mir_of_interest,
#    paste(mir_of_interest[!(mir_of_interest %in% unlist(miRNA_expression_data_dt[, gene_name]))], "-1", sep = ""),
#    paste(mir_of_interest[!(mir_of_interest %in% unlist(miRNA_expression_data_dt[, gene_name]))], "-2", sep = "")
#)


##########
## Extract the miRNA expression data related to the miRNA
##########

















############################################
## Extract the mir expression related to the mir that interact with the neurotranmission related genes
############################################

## extract the expression data related to the mir of interest
miRNA_expression_data_filtered_dt <- miRNA_expression_data_dt[gene_name %in% mir_of_interest, ]

## do the transposition of the dat table
miRNA_expression_data_filtered_dt <- transpose_datatable(
    miRNA_expression_data_filtered_dt,
    column_name = "gene_name",
    new_name = "sample_id"
)


## merge with the cluster data
miRNA_expression_data_filtered_dt <- merge(
    cluster_data_dt,
    miRNA_expression_data_filtered_dt,
    by = "sample_id",,
    all.y = T
)







## for the sample that does not have a cluster name, these samples are healthy sample so put this information into the data
miRNA_expression_data_filtered_dt[is.na(cluster), cluster := "Healthy"]


## melt the data
miRNA_expression_data_filtered_melted_dt <- melt(
    miRNA_expression_data_filtered_dt,
    id.vars = c("sample_id", "cluster"),
    measure.vars = colnames(miRNA_expression_data_filtered_dt)[!(colnames(miRNA_expression_data_filtered_dt) %in% c("sample_id", "cluster"))]
)










### do the average by cluster
#miRNA_expression_data_filtered_average_dt <- miRNA_expression_data_filtered_melted_dt[, .(average = mean(as.numeric(value), na.rm = TRUE)), by = c("cluster", "variable")]

### update the miRNA vector of interest that is contained by the expression data
#mir_of_interest <- as.vector(unique(unlist(miRNA_expression_data_filtered_average_dt[, 'variable', with = F])))

### remove the NA values from the exression profiles
#mir_of_interest <- mir_of_interest[!(is.na(mir_of_interest))]























############################################
## extract the miRNA that are significant from the deseq analysis
############################################


## set the pvalue and fold change cutoff
pval_cutoff <- as.numeric(snakemake@params[["pval_cutoff"]])
logfoldchange_cutoff <- as.numeric(snakemake@params[["logfoldchange_cutoff"]])

## extract the probes filtering with the cutoff values
DE_miRNA_results_dt <- Reduce(
    x = mir_deseq_results_dt,
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
        dt3 <- rbind(dt1, dt2)[order(log2FoldChange),]

        ## return output
        return(dt3)
    }
)


## retreive the DE mir
DE_mir_vector <- unique(unlist(DE_miRNA_results_dt[, gene_name]))

## retreive the DE mir that are significant in all analyses
DE_common_gene <- as.data.table(table(unlist(DE_miRNA_results_dt[, gene_name])))[order(N, decreasing = T),][N == 5,V1]




## retreive the average data related to the DE miRNA
#miRNA_expression_data_filtered_average_dt <- miRNA_expression_data_filtered_average_dt[variable %in% DE_mir_vector,]

#print(miRNA_expression_data_filtered_average_dt)



## filter the miRNA expression with the differentially expressed miRNA
miRNA_expression_data_filtered_melted_dt <- miRNA_expression_data_filtered_melted_dt[variable %in% DE_mir_vector,]


## filter the miRNA 
miRNA_expression_data_filtered_dt <- miRNA_expression_data_filtered_dt[
    ,
    c("sample_id", "cluster", colnames(miRNA_expression_data_filtered_dt)[colnames(miRNA_expression_data_filtered_dt) %in% DE_mir_vector]),
    with = F
]






#print(miRNA_expression_data_filtered_melted_dt[variable %in% DE_common_gene,])

#break


### write the table
#fwrite(
#    unique(miRNA_expression_data_filtered_melted_dt[order(variable, cluster),][, -c("sample_id")]),
#    paste(
#        snakemake@output[["miRNA_figure_dir"]],
#        '/',
#        "miRNA_counting",
#        '.csv',
#        sep = ""
#    ),
#    sep = ","
#)





############################################
## Filter the expression data, removing miRNA with a average lower than XXX
############################################


## calculate the average of expression per miRNA and per cluster
expression_average_per_mir <- copy(miRNA_expression_data_filtered_melted_dt[, .(expression_average_cluster = mean(as.numeric(value))), by = c('variable', 'cluster')])





## write the table
fwrite(
    expression_average_per_mir,
    paste(
        snakemake@output[["miRNA_figure_dir"]],
        '/',
        "miRNA_counting",
        '.csv',
        sep = ""
    ),
    sep = ","
)




















#miRNA_expression_data_filtered_melted_dt[, expression_average_percluster := mean(as.numeric(value)), by = c('variable', 'cluster')]


#miRNA_expression_data_filtered_melted_dt[, expression_average_permir := mean(as.numeric(value)), by = c('variable')]


#print(miRNA_expression_data_filtered_melted_dt)
#print(miRNA_expression_data_filtered_melted_dt)

#break



### calculate the average expression per miRNA
#expression_average_per_mir[, expression_average := mean(expression_average_cluster), by = c('variable')]

### retreive the miRNA that are associated with a lower expression
#expression_average_per_mir <- expression_average_per_mir[expression_average > as.numeric(snakemake@params[["expression_cutoff"]]),]
#expression_average_per_mir <- expression_average_per_mir[order(expression_average, decreasing = T),]




#DE_mir_vector <- unlist(expression_average_per_mir[,variable])

### extract the expression data associated with the lower expressed miRNA
#miRNA_expression_data_filtered_average_dt <- miRNA_expression_data_filtered_average_dt[variable %in% DE_mir_vector, ]















##########
## Set the order of the miRNA for the figures
##########


## extract the miRNA names with a specific order
DE_mir_vector2 <- unlist(expression_average_per_mir[cluster == "NT-1"][order(expression_average_cluster, decreasing = T),][,variable])

## set the order order of the variable based on the average expression value for the figures
miRNA_expression_data_filtered_melted_dt$variable <- factor(
    miRNA_expression_data_filtered_melted_dt$variable,
    levels = unique(DE_mir_vector2)
)



#print(miRNA_expression_data_filtered_melted_dt)
#break

### filter the metlted data and set the order for the figures generation
#miRNA_expression_data_filtered_melted_dt <- miRNA_expression_data_filtered_melted_dt[variable %in% DE_mir_vector2,]



## set the order order of the variable based on the average expression value for the figures
miRNA_expression_data_filtered_melted_dt$variable <- factor(
    miRNA_expression_data_filtered_melted_dt$variable,
    levels = unique(DE_mir_vector2)
)

##########
## split the data depending of the cluster rank of the miRNA expression
##########

## retreive the miRNA contained by the miRNA_expression_data_filtered_melted_dt
miRNA_vector <- unique(unlist(miRNA_expression_data_filtered_melted_dt[,variable])) 

## initialize the list that will contain all the rank associated with the miRNAs
miRNA_rank_list_overexpressed <- list()
miRNA_rank_list_underexpressed <- list()

## for each miRNA, retreive the cluster rank based on the expression
for (i_miRNA in seq(1, length(miRNA_vector))) {

    ## retreive the miRNA name
    miRNA_name <- as.character(unlist(miRNA_vector[i_miRNA]))


    ## retreive the expression data associated with the miRNA and order thed data with expression values
    miRNA_data_subset <- expression_average_per_mir[variable == miRNA_name,][order(expression_average_cluster, decreasing = T),]


    ## remove the healthy tissues
    miRNA_data_subset <- miRNA_data_subset[cluster != "Healthy",]
    

    ## retreive the cluster name associated with the higher average data
    rank1_cluster_name <- unlist(miRNA_data_subset[1, cluster])
    last_rank_cluster_name <- unlist(miRNA_data_subset[nrow(miRNA_data_subset), cluster])
    
    ## put the miRNA into the rank list
    miRNA_rank_list_overexpressed[[rank1_cluster_name]] <- append(miRNA_rank_list_overexpressed[[rank1_cluster_name]], list(miRNA_name))

    miRNA_rank_list_underexpressed[[last_rank_cluster_name]] <- append(miRNA_rank_list_underexpressed[[last_rank_cluster_name]], list(miRNA_name))

}


##########
## generate the barchart that describe the number of highest expressed genes
##########

## do the counting of mirna for each cluster
mir_counting_overexpressed <- lapply(
    miRNA_rank_list_overexpressed,
    function(x) {

        mi_RNA_number <- length(x)
        return(mi_RNA_number)
    }
)
mir_counting_underexpressed <- lapply(
    miRNA_rank_list_underexpressed,
    function(x) {

        mi_RNA_number <- length(x)
        return(mi_RNA_number)
    }
)



## create a data table containing the mir counting
mir_counting_dt_over <- data.table(
    "cluster" = names(mir_counting_overexpressed),
    "Highest-expressed miRNA" = unlist(mir_counting_overexpressed)
)
mir_counting_dt_under <- data.table(
    "cluster" = names(mir_counting_underexpressed),
    "Lowest-expressed miRNA" = unlist(mir_counting_underexpressed)
)

## merge the data together
mir_counting_dt <- merge(
    mir_counting_dt_over,
    mir_counting_dt_under,
    by = "cluster",
)[order(cluster),]

## melt data
mir_counting_dt_melted <- melt(
    mir_counting_dt,
    id.vars = "cluster",
    measure.vars = c("Highest-expressed miRNA", "Lowest-expressed miRNA")
)





## rename column
mir_counting_dt_melted <- mir_counting_dt_melted[, setnames(.SD, "variable", "miRNA Type")]

## generate the barchart
counting_barchart <- ggplot(
    data = mir_counting_dt_melted,
    aes(
        x = cluster,
        y = value,
        fill = `miRNA Type`
    )
) +
    geom_bar(stat="identity", position=position_dodge(), color = "black", size = 2) +
    scale_fill_manual(values = brewer.pal(9, "Set2")) +
    labs(
        x = "Neurotransmission Cluster",
        y = "Number of miRNAs",
        fill =''
    ) + 
    scale_y_continuous(breaks = c(seq(0, 34, 2))) + 
    theme_bw(base_size = 24)

## save the figure into a file
ggsave(
    counting_barchart,
    filename = paste(
        snakemake@output[["miRNA_figure_dir"]],
        '/',
        "miRNA_counting",
        '.svg',
        sep = ""
    ),
    device = "svg",
    height = 8,
    width = 14,
    bg = "white"
)




















## initialize the lsit of plots
plot_list <- list()

## get the cluster name into a vector
cluster_vector <- unique(unlist(miRNA_expression_data_filtered_melted_dt[,cluster]))

for (i_cluster in seq(1, length(cluster_vector))) {

    ## get the cluster name
    cluster_value <- cluster_vector[i_cluster]

    #print("==============================")
    #print(cluster_value)
    ## get the highest and lowest expressed genes based on the cluster
    lowest_expressed_genes <- unlist(miRNA_rank_list_underexpressed[[cluster_value]])
    highest_expressed_genes <- unlist(miRNA_rank_list_overexpressed[[cluster_value]])

    ## get the expression data related to this genes
    sub_expression_data_high <- expression_average_per_mir[variable %in% c(highest_expressed_genes),]
    sub_expression_data_low <- expression_average_per_mir[variable %in% c(lowest_expressed_genes),]
    sub_expression_data <- expression_average_per_mir[variable %in% c(lowest_expressed_genes, highest_expressed_genes),]

    ## set the order of the genes
    gene_order <- c(
        as.character(unique(unlist(sub_expression_data_high[cluster == cluster_value,][order(-expression_average_cluster),variable]))),
        as.character(unique(unlist(sub_expression_data_low[cluster == cluster_value,][order(expression_average_cluster),variable])))

    )

    ## generate the name of the figure
    file_name <- paste(
        "expression_curve_highest_lowest_expressed_",
        cluster_value,
        ".svg",
        sep = ""
    )

    sub_expression_data$variable <- factor(
        sub_expression_data$variable,
        levels = gene_order
    )

    ## generate the miRNA expression profile
    miRNA_expression_profiles <- ggplot(
        data = sub_expression_data,
        aes(
            x = variable,
            y = expression_average_cluster,
            color = cluster
        )
    ) +
        geom_point(alpha = 0.3) +
        geom_line(aes(group = cluster), alpha = 0.7, lwd = 2) +
        scale_color_manual(values = color_palette) +
        labs(
            x = "miRNA",
            y = "Average DESeq2 Normalized Expression"
        ) +
        theme_bw(base_size = 24) +
        theme(
            #axis.text.x = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.key.size = unit(1, 'cm')

        )


    ## add the plot into the list
    plot_list <- append(
        plot_list,
        list(miRNA_expression_profiles)
    )
    
    ggsave(
        miRNA_expression_profiles,
        filename = paste(
            snakemake@output[['miRNA_figure_dir']],
            '/',
            file_name,
            sep = ""
        ),
        device = "svg",
        height = 8,
        width = 20,
        bg = "white",
        limitsize = FALSE    
    )

}


## arrange the plots together
barchart_all <- ggarrange(
    plot_list[[3]],
    plot_list[[2]],
    plot_list[[4]],
    plot_list[[1]],
    nrow = 4,
    ncol = 1
    #heights = c(1,1),
    #widths = c(8,2)
)

## save the methylation profile average into a file
svg(
    filename = paste(
        snakemake@output[['miRNA_figure_dir']],
        '/',
        "expression_curve_highest_lowest_expressed_all",
        '.svg',
        sep = ""
    ),
    height = 25,
    width = 15
)
print(barchart_all)
dev.off()




## get the cluster name in order
cluster_vector <- sort(unique(unlist(miRNA_expression_data_filtered_melted_dt[,cluster])))


## initialize the lsit of plots
plot_list <- list()

## initialize the list that contain the gene order
gene_order_high <- c()
gene_order_low <- c()


for (i_cluster in seq(1, length(cluster_vector))) {

    ## get the cluster name
    cluster_value <- cluster_vector[i_cluster]

    print("==============================")
    print(cluster_value)
    ## get the highest and lowest expressed genes based on the cluster
    lowest_expressed_genes <- unlist(miRNA_rank_list_underexpressed[[cluster_value]])
    highest_expressed_genes <- unlist(miRNA_rank_list_overexpressed[[cluster_value]])

    ## get the expression data related to this genes
    sub_expression_data_high <- expression_average_per_mir[variable %in% c(highest_expressed_genes),]
    sub_expression_data_low <- expression_average_per_mir[variable %in% c(lowest_expressed_genes),]
    sub_expression_data <- expression_average_per_mir[variable %in% c(lowest_expressed_genes, highest_expressed_genes),]

    ## set the order of the genes
    gene_order_high <- c(
        gene_order_high,
        as.character(unique(unlist(sub_expression_data_high[cluster == cluster_value,][order(-expression_average_cluster),variable])))
    )

    gene_order_low <- c(
        gene_order_low,
        as.character(unique(unlist(sub_expression_data_low[cluster == cluster_value,][order(expression_average_cluster),variable])))

    )
}

expression_average_per_mir_high <- copy(expression_average_per_mir)
expression_average_per_mir_low <- copy(expression_average_per_mir)

expression_average_per_mir_high$variable <- factor(
    expression_average_per_mir_high$variable,
    levels = gene_order_high
)
expression_average_per_mir_low$variable <- factor(
    expression_average_per_mir_low$variable,
    levels = gene_order_low
)

## generate the miRNA expression profile
miRNA_expression_profiles_high <- ggplot(
    data = expression_average_per_mir_high,
    aes(
        x = variable,
        y = log2(expression_average_cluster),
        color = cluster
    )
) +
    geom_point(alpha = 0.3) +
    geom_line(aes(group = cluster), alpha = 0.7, lwd = 2) +
    scale_color_manual(values = color_palette) +
    labs(
        x = "miRNA",
        y = "DESeq2 Normalized Expression Average (Log2)"
    ) +
    theme_bw(base_size = 24) +
    theme(
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(1, 'cm')
    )


## generate the miRNA expression profile
miRNA_expression_profiles_low <- ggplot(
    data = expression_average_per_mir_low,
    aes(
        x = variable,
        y = log2(expression_average_cluster),
        color = cluster
    )
) +
    geom_point(alpha = 0.3) +
    geom_line(aes(group = cluster), alpha = 0.7, lwd = 2) +
    scale_color_manual(values = color_palette) +
    labs(
        x = "miRNA",
        y = "Average DESeq2 Normalized Expression"
    ) +
    theme_bw(base_size = 24) +
    theme(
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(1, 'cm')

    )



## arrange the plots together
barchart_all <- ggarrange(
    miRNA_expression_profiles_high,
    miRNA_expression_profiles_low,
    nrow = 2,
    ncol = 1
    #heights = c(1,1),
    #widths = c(8,2)
)

## save the methylation profile average into a file
svg(
    filename = paste(
        snakemake@output[['miRNA_figure_dir']],
        '/',
        "expression_curve_highest_lowest_expressed_allmir",
        '.svg',
        sep = ""
    ),
    height = 14,
    width = 20
)
print(barchart_all)
dev.off()




############################################
## Generate the expression barplots of expression average for all the DE miRNA
############################################


## extract the DE miRNAs
#mir_vector <- c(
#    unlist(miRNA_rank_list_overexpressed, recursive = T),
#    unlist(miRNA_rank_list_underexpressed, recursive = T)
#)

#print(mir_vector)
#break
## extract the expression data related to the mir of interest
#mir_expression_data_dt <- miRNA_expression_data_filtered_melted_dt[variable %in% mir_vector,]

## get the expression for each miRNA of interest and for each sample
mir_expression_of_interest <- copy(miRNA_expression_data_filtered_melted_dt)

## calculate the mean for each miRNA
mir_expression_of_interest[, expression_average := mean(as.numeric(value)), by = variable] 
mir_expression_of_interest[, standard_deviation := sd(as.numeric(value)), by = variable] 

## remove the sample id and the expression values and keep the average expression and the standard deviation
mir_expression_of_interest <- unique(mir_expression_of_interest[, !c("sample_id", "value"), with = F])

##########
## generate the figure
##########

## order by the mean
mir_expression_of_interest$variable <- factor(
    mir_expression_of_interest$variable,
    levels = unique(unlist(mir_expression_of_interest[order(-expression_average),variable]))
)

## generate the figures
average_expression_barchart <- ggplot(
    data = mir_expression_of_interest,
    aes(
        x = variable,
        y = expression_average
    )
) + 
    geom_bar(stat="identity", position=position_dodge(), color = "black", fill = "grey", size = 2) +
    geom_errorbar(aes(ymin=expression_average-standard_deviation, ymax=expression_average+standard_deviation), width=.2, position=position_dodge(.9), size = 1.5) + 
    geom_hline(yintercept = mean(mir_expression_of_interest[, expression_average]), color = "red", size = 3) +
    #scale_fill_manual(values = '#000000') +
    labs(
        x = "",
        y = "Average DESeq2 Normalized Expression"
    ) +
    theme_bw(base_size = 60) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none",
        legend.key.size = unit(1.5, 'cm')
    )

## save the figure into a file
ggsave(
    average_expression_barchart,
    filename = paste(
        snakemake@output[["miRNA_figure_dir"]],
        #'data',
        '/',
        "miRNA_expression_average_bartchart",
        '.svg',
        sep = ""
    ),
    device = "svg",
    height = 25,
    width = 60,
    bg = "white",
    limitsize = FALSE  
)





































############################################
## generate a barchart with the best miRNA whose expression are higher than the average expression
############################################


## copy the expression data related to all the mir
mir_expression_of_interest_cluster <- copy(miRNA_expression_data_filtered_melted_dt)

## do the expression average for each cluster
mir_expression_of_interest_cluster[, cluster_expression_average := mean(as.numeric(value)), by = c("variable", "cluster")]
mir_expression_of_interest_cluster[, standard_deviation := sd(as.numeric(value)), by = c("variable", "cluster")]


## get the miRNA of interest whose expression are higher than the average 
high_expressed_mir_vector <- unique(unlist(mir_expression_of_interest[expression_average > mean(mir_expression_of_interest[, expression_average]) , variable]))



## get the expression related to the high expressed miRNA
high_expressed_mir_per_cluster <- mir_expression_of_interest_cluster[variable %in% high_expressed_mir_vector,]
high_expressed_mir_per_cluster <- unique(high_expressed_mir_per_cluster[, !c("sample_id", "value")], with = F)[order(variable), ]


print(mean(mir_expression_of_interest[, expression_average]))
print(unique(unlist(high_expressed_mir_per_cluster[,variable])))



## get the top 3 of the miRNA
top_mir <- mir_expression_of_interest_cluster[, mean(as.numeric(value)), by = "variable"][order(V1, decreasing = T),][1:4, variable]

high_expressed_mir_per_cluster <- mir_expression_of_interest_cluster[variable %in% top_mir,]

##########
## Generate the figures
##########

## order the cluster categorie for the vizualization
high_expressed_mir_per_cluster$cluster <- factor(
    high_expressed_mir_per_cluster$cluster,
    levels = c(
        "NT-1",
        "NT-2",
        "NT-3",
        "NT-4",
        "Healthy"
    )
)

## generate the figures
high_expressed_barchart <- ggplot(
    data = high_expressed_mir_per_cluster, 
    aes(
        x = variable,
        y = (cluster_expression_average),
        fill = cluster
    )
) +
    geom_bar(
        stat="identity",
        position=position_dodge(),
        color = "black",
        size = 2
    ) +
    geom_errorbar(aes(ymin=cluster_expression_average-standard_deviation, ymax=cluster_expression_average+standard_deviation), width=.2, position=position_dodge(.9), size = 1.5) + 
    scale_fill_manual(values = brewer.pal(9, "Set1")) +
    labs(
        x = "Top 4 High Expressed miRNAs",
        y = "DESeq2 Normalized Expression Average"
    ) + 
    ylim(0,25) +
    theme_bw(base_size = 24) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        #legend.position = "none",
        legend.key.size = unit(1.5, 'cm')
    )

## save the figure into a file
ggsave(
    high_expressed_barchart,
    filename = paste(
        snakemake@output[["miRNA_figure_dir"]],
        '/',
        "high_expressedmiRNA_expression_average_bartchart",
        '.svg',
        sep = ""
    ),
    device = "svg",
    height = 10,
    width = 15,
    bg = "white"
)







































































###############################
## TO DELETE
###############################









### calculate the average for all the miRNAs of interest
#mir_expression_data_dt <- mir_expression_data_dt[, .(value = mean(as.numeric(value), na.rm = T)), by = "variable"][order(value, decreasing = T),]

### do the mean of the expression
#mean_expression <- mean(as.numeric(unlist(mir_expression_data_dt[,value])))
#print(mir_expression_data_dt)
#cat("The expression average for the DE miRNA:", mean_expression, "\n")


### set the order
#mir_expression_data_dt$variable <- factor(
#    mir_expression_data_dt$variable,
#    levels = unlist(mir_expression_data_dt$variable)
#)

### generate the figure
#average_expression_barchart <- ggplot(
#    data = mir_expression_data_dt,
#    aes(
#        x = variable,
#        y = value
#    )
#) + 
#    geom_bar(stat="identity", position=position_dodge(), color = "black", fill = "grey", size = 2) +
#    geom_hline(yintercept = mean_expression, color = "red", size = 3) +
#    #scale_fill_manual(values = '#000000') +
#    labs(
#        x = "",
#        y = "Average Expression levels (RPKM)"
#    ) +
#    theme_bw(base_size = 60) +
#    theme(
#        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#        legend.position = "none",
#        legend.key.size = unit(1.5, 'cm')
#    )

### save the figure into a file
#ggsave(
#    average_expression_barchart,
#    filename = paste(
#        snakemake@output[["miRNA_figure_dir"]],
#        #'data',
#        '/',
#        "miRNA_expression_average_bartchart",
#        '.svg',
#        sep = ""
#    ),
#    device = "svg",
#    height = 25,
#    width = 60,
#    bg = "white",
#    limitsize = FALSE  
#)



###########
### Generate the expression barplots for each rank cluster name
###########

### initialize the list that will contain all the plots
#plot_list <- list()

#for (i_miRNA_rank_list in seq(1, length(miRNA_rank_list_overexpressed))) {

#    ## retreive the cluster name first ranked
#    cluster_name_rank1 <- names(miRNA_rank_list_overexpressed[i_miRNA_rank_list])

#    ## retreive the miRNA associated with the cluster name rank1
#    miRNA_rank1_vector <- unlist(miRNA_rank_list_overexpressed[[cluster_name_rank1]], recursive = T)

#    ## extract the expression data associated with the miRNA ranked 1
#    miRNA_data_subset_expression <- expression_average_per_mir[variable %in% miRNA_rank1_vector, ]

#    ## extract the miRNA ordered with the expression data based on the cluster name rank 1
#    miRNA_vector_ordered <- expression_average_per_mir[cluster == cluster_name_rank1][order(expression_average_cluster, decreasing = T),][, variable]

#    ## set the order order of the variable based on the expression_average_cluster expression expression_average_cluster for the figures
#    miRNA_data_subset_expression$variable <- factor(
#        miRNA_data_subset_expression$variable,
#        levels = unique(miRNA_vector_ordered)
#    )

#    ## set the order of the cluster, putting the Healthy condition at the end
#    cluster_name_vector <- unique(unlist(miRNA_data_subset_expression$cluster))
#    cluster_name_vector <- sort(cluster_name_vector[cluster_name_vector != "Healthy"])
#    cluster_name_vector <- c(cluster_name_vector, "Healthy")

#    ## set the order of the cluster
#    miRNA_data_subset_expression$cluster <- factor(
#        miRNA_data_subset_expression$cluster,
#        levels = cluster_name_vector
#    )

#    ## generate hte barchart for each miRNA
#    miRNA_expression_bartchart <- ggplot(
#        data = miRNA_data_subset_expression,
#        aes(
#            x = variable,
#            y = log2(expression_average_cluster),
#            fill = cluster
#        )
#    ) +
#        geom_bar(stat="identity", position=position_dodge(), color = "black", size = 1) +
#        #geom_boxplot() +
#        #geom_violin(size = 2, color = "black") +
#        #scale_y_continuous(limits=c(0, 300)) +
#        scale_fill_manual(values = color_palette) +
#        labs(
#            x = "",
#            y = "Log2 Average Expression Levels (RPKM)"
#        ) +
#        theme_bw(base_size = 70) +
#        theme(
#            #axis.text.x = element_blank(),
#            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#            legend.position = "none",
#            legend.key.size = unit(1.5, 'cm')

#        )

#    ## add the plot into the list
#    plot_list <- append(
#        plot_list,
#        list(miRNA_expression_bartchart)
#    )


#    ## save the figure into a file
#    ggsave(
#        miRNA_expression_bartchart,
#        filename = paste(
#            snakemake@output[["miRNA_figure_dir"]],
#            #'data',
#            '/',
#            "miRNA_expression_bartchart_",
#            cluster_name_rank1,
#            '.svg',
#            sep = ""
#        ),
#        device = "svg",
#        height = 25,
#        width = 50,
#        bg = "white",
#        limitsize = FALSE  
#    )

#}



### arrange the plots together
#barchart_all <- ggarrange(
#    plot_list[[3]],
#    plot_list[[2]],
#    plot_list[[1]],
#    plot_list[[4]],
#    nrow = 4,
#    ncol = 1
#    #heights = c(1,1),
#    #widths = c(8,2)
#)


### arrange the plots together
##barchart_all <- ggarrange(
##    plot_list[[2]],
##    plot_list[[3]],
##    plot_list[[1]],
##    plot_list[[4]],
##    nrow = 2,
##    ncol = 2,
##    heights = c(1,1),
##    widths = c(8,2)
##)

### arrange to merge plots togeter
##barchart_all <- ggarrange(
##    barchart_all,
##    average_expression_barchart,
##    nrow = 2
##)

### save the methylation profile average into a file
#svg(
#    filename = paste(
#        snakemake@output[['miRNA_figure_dir']],
#        '/',
#        "miRNA_expression_bartchart_all",
#        '.svg',
#        sep = ""
#    ),
#    height = 60,
#    width = 30
#)
#print(barchart_all)
#dev.off()




## generate the miRNA expression barchart for all miRNA together
#miRNA_expression_profiles <- ggplot(
#    data = miRNA_expression_data_filtered_average_dt,
#    aes(
#        x = variable,
#        y = average,
#        fill = cluster
#    )
#) +
#    geom_bar(stat="identity", position=position_dodge()) +
#    #geom_boxplot() +
#    #scale_y_continuous(limits=c(0, 100)) +
#    scale_fill_manual(values = color_palette) +
#    labs(
#        x = "miRNA",
#        y = "Average Methylation Levels\n(Beta Values)"
#    ) +
#    theme_bw(base_size = 60) +
#    theme(
#        #axis.text.x = element_blank(),
#        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#        legend.key.size = unit(1.5, 'cm')

#    )

#ggsave(
#    miRNA_expression_profiles,
#    filename = paste(
#        #snakemake@output[['methylation_figure_dir']],
#        'data',
#        '/',
#        "miRNA_expression_profiles",
#        '.svg',
#        sep = ""
#    ),
#    device = "svg",
#    height = 30,
#    width = 75,
#    bg = "white",
#    limitsize = FALSE  
#)



















##########
## Generate the miRNA expression profiles
##########


#print(expression_average_per_mir)


test <- copy(expression_average_per_mir)
print(test)

## dcast the data
test <- dcast(
    test,
    variable ~ cluster,
    value.var = "expression_average_cluster"
)

print(test)
setorderv(test, c("NT-2", "NT-1"),  order = c(-1, 1))
gene_order <- unlist(test[, variable])




expression_average_per_mir$variable <- factor(
    expression_average_per_mir$variable,
    levels = unique(unlist(expression_average_per_mir[cluster == "NT-1",][order(-expression_average_cluster),variable])),
    #levels = gene_order
)



## generate the miRNA expression profile
miRNA_expression_profiles <- ggplot(
    data = expression_average_per_mir,
    aes(
        x = variable,
        y = expression_average_cluster,
        color = cluster
    )
) +
    geom_point(alpha = 0.3) +
    geom_line(aes(group = cluster), alpha = 0.7, lwd = 5) +
    #geom_bar() +
    scale_color_manual(values = color_palette) +
    labs(
        x = "miRNA",
        y = "Average DESeq2 Normalized Expression"
    ) +
    theme_bw(base_size = 60) +
    theme(
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(1.5, 'cm')

    )

ggsave(
    miRNA_expression_profiles,
    filename = paste(
        snakemake@output[['miRNA_figure_dir']],
        #'data',
        '/',
        "miRNA_expression_profiles",
        '.svg',
        sep = ""
    ),
    device = "svg",
    height = 30,
    width = 50,
    bg = "white",
    limitsize = FALSE    
)



##########
## Generate the average miRNA expression boxplot
##########


## generate the boxplot of the miRNA expression per groups
miRNA_expression_boxplot <- ggplot(
    data = miRNA_expression_data_filtered_melted_dt,
    aes(
        x = cluster,
        y = log2(as.numeric(value)),
        fill = cluster
    )
) +
    #geom_boxplot(size = 2, color = 'black', notch = T) +
    geom_violin(size = 2, color = "black") +
    scale_fill_manual(values = color_palette) +
    labs(
        x = "Cluster",
        y = "miRNA Expression Level"
    ) +
    theme_minimal(base_size = 24)


ggsave(
    miRNA_expression_boxplot,
    filename = paste(
        snakemake@output[['miRNA_figure_dir']],
        #'data',
        '/',
        "miRNA_expression_boxplot",
        '.svg',
        sep = ""
    ),
    device = "svg",
    height = 8,
    width = 10,
    bg = "white"    
)


##########
## Generate a heatmap based on the miRNA expression
##########

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
create_heatmap_miRNA_expression <- function(
    miRNA_expression_matrix,
    clinical,
    variables,
    path_file,
    clustering_distance_rows = "minkowski",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "pearson",
    clustering_method_columns = "ward.D2"
) {

    ## generate the heatmap file
    generation_heatmap <- FALSE

    ## do the transpose of the matrix
    miRNA_expression_matrix <- t(miRNA_expression_matrix)
    
    if (!(missing(path_file))) {
        generation_heatmap <- TRUE
    }

    ## extract and make the annotations matrix
    annotation <- create_annotation(variables, clinical)

    ## generate the heatmap
    ht <- Heatmap(
        miRNA_expression_matrix,
        border = TRUE,
        #column_title = column_title,
        #row_title = row_title,
        #column_title_gp = column_title_gp,
        #row_title_gp = row_title_gp,
        show_column_names = FALSE,
        show_row_names = FALSE,


        #cluster_columns = FALSE,

        
        ## legends for the heatmap
        top_annotation = annotation,
        heatmap_legend_param = list(
            title = "miRNA Expression Level",
            title_position = "leftcenter-rot",
            #at = c(0, 0.5, 1),
            grid_width = unit(0.75, "cm"),
            legend_height = unit(6, "cm"),
            border = "black"
        ),

        ## number of clusters
        column_split = 5,
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



#####
## create the heatmap for each condition
######


## get the cluster of interest
cluster_of_interest <- unique(names(miRNA_rank_list_underexpressed))
print(cluster_of_interest)

## for each cluster of interest, extracting the data of interest
for (i_cluster in seq(1, length(cluster_of_interest))) {

    ## get the cluster of interest
    cluster_value <- cluster_of_interest[i_cluster]

    ## copy the mirna of the cluster of interest
    miRNA_expression_data_subset <- copy(miRNA_expression_data_filtered_dt[cluster == cluster_value,])

    ## transform the mir expression data to matrix
    miRNA_expression_matrix <- as.matrix(miRNA_expression_data_subset[, !c("cluster", "sample_id"), with = F])

    ## transform all the values to numeric values
    miRNA_expression_matrix <- apply(miRNA_expression_matrix, c(1,2), function(x) log2(as.numeric(x)))

    ## extract the information cluster associated with the samples
    cluster_data <- miRNA_expression_data_subset[, c("cluster", "sample_id"), with = F]

    ## call the function for the heatmap generation
    create_heatmap_miRNA_expression(
        miRNA_expression_matrix,
        cluster_data,
        c("cluster"),
        clustering_distance_columns = "pearson",
        path_file = paste(
            snakemake@output[['miRNA_figure_dir']],
            #"data",
            '/',
            "miRNA_expression_heatmap_",
            cluster_value,
            '.svg',
            sep = ""
        ),
    )

    next
    


    ### get the genes that are overexpressed or underepressed in this cluster
    #overexpressed_genes <- unlist(miRNA_rank_list_overexpressed[[cluster_value]])
    #underexpressed_genes <- unlist(miRNA_rank_list_underexpressed[[cluster_value]])

    ### extracting the mir expression of interest
    #expression_subset_over <- miRNA_expression_data_filtered_dt[
    #    #cluster == cluster_value,
    #    ,
    #    c(
    #        "sample_id",
    #        "cluster",
    #        overexpressed_genes
    #    ),
    #    with = F
    #]

    #expression_subset_under <- miRNA_expression_data_filtered_dt[
    #    #cluster == cluster_value,
    #    ,
    #    c(
    #        "sample_id",
    #        "cluster",
    #        underexpressed_genes
    #    ),
    #    with = F
    #]
    
    ######
    ### Formating the data for generating the heatmap 
    ######
    

    ### transform the mir expression data to matrix
    #expression_subset_matrix_over <- as.matrix(expression_subset_over[, !c("cluster", "sample_id"), with = F])
    #expression_subset_matrix_under <- as.matrix(expression_subset_under[, !c("cluster", "sample_id"), with = F])

    ### transform all the values to numeric values
    #expression_subset_matrix_over <- apply(expression_subset_matrix_over, c(1,2), as.numeric)
    #expression_subset_matrix_under <- apply(expression_subset_matrix_under, c(1,2), as.numeric)

    ### extract the information cluster associated with the samples
    #cluster_data_over <- expression_subset_over[, c("cluster", "sample_id"), with = F]
    #cluster_data_under <- expression_subset_under[, c("cluster", "sample_id"), with = F]

    ### put the heatnmap into the list
    #heatmap_list_over <- create_heatmap_miRNA_expression(
    #    expression_subset_matrix_over,
    #    cluster_data_over,
    #    c("cluster"),
    #    path_file = paste(
    #        snakemake@output[['miRNA_figure_dir']],
    #        #"data",
    #        '/',
    #        "miRNA_expression_heatmap_over",
    #        '.svg',
    #        sep = ""
    #    ),
    #)

    #heatmap_list_under <- create_heatmap_miRNA_expression(
    #    expression_subset_matrix_under,
    #    cluster_data_under,
    #    c("cluster"),
    #    path_file = paste(
    #        snakemake@output[['miRNA_figure_dir']],
    #        #"data",
    #        '/',
    #        "miRNA_expression_heatmap_under",
    #        '.svg',
    #        sep = ""
    #    ),
    #)

    #print("okkk")
    #Sys.sleep(10000)


    #print("---------------------------------")
    ##print(expression_subset_over[1:4,1:4])
    #print(dim(expression_subset_over))
    ##print(cluster_value)
    ##print(overexpressed_genes)
    ##print(underexpressed_genes)



}

print("okkkkk")
Sys.sleep(100000)
break
































print(miRNA_expression_data_filtered_dt[1:4,1:4])

miRNA_expression_data_filtered_dt <- miRNA_expression_data_filtered_dt[order(cluster),]



## transform the mir expression data to matrix
miRNA_expression_matrix <- as.matrix(miRNA_expression_data_filtered_dt[, !c("cluster", "sample_id"), with = F])


## transform all the values to numeric values
#miRNA_expression_matrix <- apply(miRNA_expression_matrix, c(1,2), function(x) log2(as.numeric(x)))
miRNA_expression_matrix <- apply(miRNA_expression_matrix, c(1,2), function(x) (as.numeric(x)))

## extract the information cluster associated with the samples
cluster_data <- miRNA_expression_data_filtered_dt[, c("cluster", "sample_id"), with = F]






## call the function for 
create_heatmap_miRNA_expression(
    miRNA_expression_matrix,
    cluster_data,
    c("cluster"),
    clustering_distance_columns = "pearson",
    path_file = paste(
        snakemake@output[['miRNA_figure_dir']],
        #"data",
        '/',
        "miRNA_expression_heatmap",
        '.svg',
        sep = ""
    ),
)



print("finiiiiiiiiiiiiiiiiiiiiiii")
Sys.sleep(100000)