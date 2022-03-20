
###################################################
###################################################

##### ANALYSE THE GENES WHICH ARE DIFFERENTIAL EXPRESSED BETWEEN THE CLUSTERS

###################################################
###################################################


##### LOAD THE LIBRARIES
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(UpSetR)
library(EnhancedVolcano)
library(egg)
library(ggpubr)

###################################################
###################################################


##### LOAD THE DATA

## source the r utils
source(snakemake@params[["utils_R"]])

## load the wildcard
wcards <- snakemake@params[["wcards"]]

## load the FPKM data that are associated with the genes of interest
FPKM_data_dt <- fread(
    "data/expression_data/TCGA_expression_data_fpkm_cancer",
    sep = ","
)

## load the ensembl id data table
id_data <- fread(
    snakemake@params[["ensembl_id"]],
    sep = ","
)

## load the ensembl id
ensembl_id <- fread(
    input = snakemake@params[["ensembl_id"]],
    sep = ","
)

## read the healthy fpkm data
FPKM_data_healthy_dt <- fread(
    snakemake@input[["normal_fpkm_expression_data"]],
    sep = ","
)

## transform the ensembl id to genes
FPKM_data_healthy_dt[, genes := ensemblID_to_geneSymbole(genes, ensembl_id)]


## load the cluster data
clusters <- fread(
    snakemake@input[["cluster_group"]],
    sep = ","
)

## retreive the cluster names of the data
cluster_names <- unique(unlist(clusters[, cluster]))

expression_data <- fread(
    snakemake@input[["normalized_expression_data_expressionFiltered"]],
    sep = ","
)

## load the directory
DEgenes_directory <- snakemake@input[["DEgene_results_dir"]]

## retreive the file names
file_names <- list.files(DEgenes_directory)

## retreive only the file names that describe just the filtered results
file_names <- grep(file_names, pattern = "filtered", value = TRUE)

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

## rename the names of each data removing the filtered_ word from the name
names(merged_file) <- str_replace_all(file_names, pattern = "filtered_", replacement = "")

## rename the names of each sublist of data
names(merged_file) <- str_match(
    names(merged_file),
    pattern = paste(
        wcards,
        "_(.*?).csv",
        sep = ""
    )
)[,2]

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

## Set the pvalue cutoff and the fold change cutoff
pval_cutoff <- 0.001
FCcutoff <- 1

#####################################################################################
## TODO : do the intersection of each group independent to the number of pathways in input

# ## do the combination of each data
# pathway_combinations <- combn(
#     x = names(KEGG_genes_list),
#     m = 2 
# )

# for (i in seq(1, ncol(pathway_combinations), 1)) {

#     ## extract the i-th combination
#     combination <- pathway_combinations[, i]
# }
#####################################################################################

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

## create the directory of volcanoplot outputs
dir.create(snakemake@output[["volcano_plot_directory"]])


###################################################
###################################################


##### GENERATE THE VOLCANOPLOTS FOR EACH RESULTS OF DESEQ2

## for each data in the merged_file
for (i in seq(1, length(merged_file), 1)) {

    # print("###############################")

    ## extract the subset_data of the merged_file list
    subset_data <- as.data.table(merged_file[[i]])

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

    ## intialize the cutoff of the genes

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
    theme_minimal(base_size = 34) +
    labs(color='Pathways') +
    labs(subtitle = title)

    ## generate the file path for saving
    file_path <- paste(
        snakemake@output[["volcano_plot_directory"]],
        "/volcano_plot_",
        subdata_name,
        ".png",
        sep = ""
    )

    ## save the volcano plot in the diretory
    # svg(
    #     filename = file_path,
    #     width = 14,
    #     height = 8
    # )
    png(
        filename = file_path,
        width = 1000,
        height = 600
    )
    print(volca_plot)
    dev.off()
    

}

# Sys.sleep(10000)
###################################################
###################################################


##### GENERATE UPSETPLOTS


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


## get the unique differential expressed genes
all_significant_gene <- unlist(significant_gene_list, recursive = T)
print(length(all_significant_gene))


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

print(upset_plot)

## save the uptset plot into the directory
svg(
    filename = snakemake@output[["upset_plot_DEgenes"]],
    width = 16,
    height = 7,
)
print(upset_plot)
dev.off()


###################################################
## FUNCTIONS
###################################################


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

    ## update the order
    expression_data_melted_average$variable <- factor(
        expression_data_melted_average$variable,
        levels = gene_order
    )

    ## return the formated expression data 
    return(expression_data_melted_average)
}

#####
## Function for the generation of the expression profil
#####

## function for generate the expression profile plot
expression_profile_plot <- function(
    expression_data_input,
    cluster_name_vector_input,
    color_table_input
) {

    ## generate the expression profile plot
    p <- ggplot(
        data = expression_data_input[ cluster %in% cluster_name_vector_input,],
        aes(
            x = variable,
            y = value,
            color = cluster
        )
    ) + 
    theme_minimal(base_size = 90) + 
    geom_point(alpha = 0) +
    geom_line(aes(group = sample_id), alpha = 0.1, lwd = 1) +

    geom_point(data = expression_data_input[cluster == "mean_profile",],  aes(x = variable, y = value, color = cluster), color = "black", size = 5, alpha = 0) +
    geom_line(data = expression_data_input[cluster == "mean_profile",],  aes(x = variable, y = value, color = cluster, group = sample_id), lwd = 3, color = "black") +
    theme_classic(base_size = 50) +
    labs(
        x = "Genes",
        y = "Normalized Expression"
    ) +
    scale_y_continuous(limits=c(0, 2.5)) + # set the scale of the y axe
    theme(
        #axis.text.x=element_blank(),
        legend.key.size = unit(1.5, 'cm')
    ) +
    scale_color_manual(
        values = color_table_input[cluster %in% cluster_name_vector_input,color] # set the color
    )

    ## generate the file_path
    file_path <- paste(
        snakemake@output[['expression_profile_plot_directory']],
        "/expression_profil_",
        paste(cluster_name_vector_input, collapse = "_"), # if the vector contain more than two values
        ".jpeg",
        sep = ""

    )
    
    ## save the plot
    ggsave(
        plot = p,
        filename = file_path,
        device = "jpeg",
        height = 15,
        width = 50,
        limitsize = F
    )

    print("Plot saved!")

    return(p)
}

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
            theme_bw(base_size = 50) +
            labs(
                subtitle = title_name,
                x = "Genes",
                y = "Normalized Expression"
            ) +
            # scale_y_continuous(limits=c(2.5, 18.5)) + # set the scale of the y axe
            #scale_y_continuous(limits = c(0.50, 2)) + # set the scale of the y axe
            theme(
                axis.text.x=element_blank(),
                #axis.line.x=element_line(),
                #axis.line.y=element_line(),
                legend.key.size = unit(1.5, "cm")
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
        height = 15,
        width = 75,
        limitsize = F
    )

    print("Plot saved!")

    return(p)
}


###################################################
## Create the directory that will contain the expression profil plots
###################################################


## create the directory that will contain all the expression profile figures
dir.create(snakemake@output[['expression_profile_plot_directory']])


###################################################
## Extract the genes of interest
###################################################


## take the unique genes of all the significant gene lists
significant_gene_vector <- unique(unlist(significant_gene_list, recursive = FALSE))


###################################################
## for eachmetabolic pathway, generate the expression profil plots
###################################################

## initialize the list that will contain the plots
expression_plot_list <- list()

## for each metabolic pathway gene set contained in the KEGG gene list
for (i_gene_set in seq(1, length(KEGG_genes_list))) {
    
    ## retreive the metabolic names
    kegg_name <- names(KEGG_genes_list[i_gene_set])
    
    ## retreive the genes related to the kegg pathway
    kegg_gene_set <- unlist(KEGG_genes_list[i_gene_set], use.names = FALSE)

    ## do the intersection with the significant genes that we identified
    kegg_gene_set_significant <- intersect(
        significant_gene_vector,
        kegg_gene_set
    )
    
    ## for each significant kegg genes, extract the expression data related to these genes
    kegg_gene_set_expression <- get_expression_data(
        kegg_gene_set_significant,
        expression_data
    )

    ## set colors that will be used for the generation of the plots
    color_palette <- brewer.pal(9, "Set1")

    ## name of the clusters that we want to do the expression profile plot
    cluster_name_vector <- (names(table(kegg_gene_set_expression[, "cluster"])))

    ## remove the mean profile
    cluster_name_vector <- cluster_name_vector[cluster_name_vector!= "mean_profile"]

    ## create a datatable containing the association of colors for each cluster name
    color_data <- data.table(
        cluster = cluster_name_vector,
        color = color_palette[1:(length(cluster_name_vector))]
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
        snakemake@output[['expression_profile_plot_directory']],
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
        snakemake@output[['expression_profile_plot_directory']],
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
## Generate the expression plot for the significant gene that were retreived
###################################################

## extract the expression data from the genes of interest
expression_data_melted_average <- get_expression_data(
    significant_gene_vector,
    expression_data
)

#####
## Initialize the color for the figures
#####

## colors that will be used for the generation of the plots
color_palette <- brewer.pal(9, "Set1")

## name of the clusters that we want to do the expression profile plot
cluster_name_vector <- (names(table(expression_data_melted_average[, "cluster"])))

## create a datatable containing the association of colors for each cluster name
color_data <- data.table(
    cluster = cluster_name_vector,
    color = color_palette[1:(length(cluster_name_vector))]
)

#####
## Generate the average expression plot
#####

## create the average expression profile plot
plot_all_expression_profile_average <- expression_profile_plot_average(
    expression_data_input = expression_data_melted_average,
    cluster_name_vector_input = cluster_names,
    color_table_input = color_data
)

## save the expression profile average into a file
svg(
    filename = paste(
        snakemake@output[['expression_profile_plot_directory']],
        '/expression_profil_all_cluster.svg',
        sep = ''
    ),
    height = 20,
    width = 40
)
print(plot_all_expression_profile_average)

dev.off()


######
### Extract the genes that are in common with all the differential expression analysis
######

## take the significant genes in common in all the significant gene lists
significant_gene_vector_in_common <- Reduce(intersect, significant_gene_list)

## extract the expression data from the genes of interest
expression_data_melted_average <- get_expression_data(
    significant_gene_vector_in_common,
    expression_data
)

## create the average expression profile plot
plot_all_expression_profile_average <- expression_profile_plot_average(
    expression_data_input = expression_data_melted_average,
    cluster_name_vector_input = cluster_names,
    color_table_input = color_data
)

## save the expression profile average into a file
svg(
    filename = paste(
        snakemake@output[['expression_profile_plot_directory']],
        '/expression_profil_all_cluster_genes_in_commun.svg',
        sep = ''
    ),
    height = 10,
    width = 30
)
print(plot_all_expression_profile_average)
dev.off()


#####
## Generate box plot of the
#####


## transform the ensembl id to gene_name
gene_vector <- unlist(FPKM_data_dt[,genes])

## remove the number after the "." from the ensembl id
gene_vector <- str_split(
    gene_vector,
    pattern = "[.]"
)
gene_vector <- sapply(gene_vector, function(x) return(x[[1]][1]))

## transform ensemb id to gene name 
gene_vector <- map_id(
    gene_vector,
    id_data,
    'gene_id',
    'gene_name'
)

## replace the ensemhl id to gene name
FPKM_data_dt[, gene_name := gene_vector]

## remove the column of the ensembl id
FPKM_data_dt[, genes := NULL]

## extract the FPKM data related to the signifcant genes
FPKM_data_dt <- FPKM_data_dt[gene_name %in% significant_gene_vector_in_common, ] 

## do the transposition
FPKM_data_dt <- transpose_datatable(
    FPKM_data_dt,
    column_name = "gene_name",
    new_name = "sample_id"
)

print(FPKM_data_dt)

## merge the FPKM data with the cluster data
FPKM_data_dt <- merge(
    FPKM_data_dt,
    clusters,
    by = 'sample_id'
)


print(significant_gene_vector_in_common)
print(FPKM_data_dt)



#########
## extract the healthy data
#########

## extract the expression data of the genes of interest from
FPKM_data_healthy_dt <- FPKM_data_healthy_dt[genes %in% significant_gene_vector_in_common,]

print(FPKM_data_healthy_dt)

## do the transposition
FPKM_data_healthy_dt <- transpose_datatable(
    FPKM_data_healthy_dt,
    column_name = "genes",
    new_name = "sample_id"
)

## Add the cluster name
FPKM_data_healthy_dt[, cluster := "Healthy"]


## merging the healthy controle data with the cancer data
FPKM_data_dt <- rbind(
    FPKM_data_dt,
    FPKM_data_healthy_dt
)


##########
## Do the statistical analysis for expression
##########

## do the statistical analysis for the three genes comparing each cluster
statistical_results <- lapply(
    significant_gene_vector_in_common,
    function(x) statistical_analysis_quantitative(
        FPKM_data_dt,
        x,
        "cluster"
    )
)

## merge the statistical sub datatable into a same data table
statistical_results <- Reduce(
    x = statistical_results,
    f = function(x, y) rbind(x, y)
)

## write the data table that contain the statistical results
fwrite(
    statistical_results,
    paste(
        snakemake@output[['expression_profile_plot_directory']],
        '/statistical_analysis_common_DE_genes.csv',
        sep = ''
    ),
    sep = ","
)


## melt the data
FPKM_data_dt <- melt(
    FPKM_data_dt,
    id.vars = c('sample_id', 'cluster'),
    measure.vars = significant_gene_vector_in_common
)

## order the cluster name
FPKM_data_dt$cluster <- factor(
    FPKM_data_dt$cluster,
    levels = c(
        "NT-1",
        "NT-2",
        "NT-3",
        "NT-4",
        "Healthy"
    )
)

## set the color
colors <- brewer.pal(9, "Set1")

## generate the boxplot
boxplot_expression <- ggplot(
    data = FPKM_data_dt,
    aes(
        x = variable,
        y = log2(as.numeric(value)),
        fill = cluster
    )
) +
    geom_boxplot(size = 2, notch = T) +
    labs(
        x = "Neurotransmission-Related Genes",
        y = "Log2(FPKM)"
    ) +
    scale_fill_manual(values = colors) +
    theme_bw(base_size = 24) +
    theme(
        legend.key.size = unit(1, "cm")
    )


## save the plot
ggsave(
    plot = boxplot_expression,
    filename = paste(
    snakemake@output[['expression_profile_plot_directory']],
    '/boxplot_expression_profile_commun_DE_genes.svg',
    sep = ''
    ),
    device = "svg",
    height = 8,
    width = 14
)



###############################################
## Write the data into files
###############################################


#####
## write the name of the genes that were used for the profil generation
#####

## write the genes that are differentially expressed
cat(
    significant_gene_vector,
    file = snakemake@output[['DE_genes']],
    sep = "\n"
)

## write the genes that are differentially expressed and in common between all the DE analyis
cat(
    significant_gene_vector_in_common,
    file = snakemake@output[['DE_genes_in_common']],
    sep = "\n"
)




