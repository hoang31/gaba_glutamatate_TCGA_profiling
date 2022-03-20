

###################################################
## Validation of results using CGGA database
###################################################


###################################################
## Load libraries
###################################################


library("RColorBrewer")
library("ComplexHeatmap")
suppressMessages(library("tidyverse"))
library("survival")
library("survminer")
library("data.table")


###################################################
## Load the data
###################################################


## load the expression data
expression_data_CGGA <- fread(
    file = snakemake@input[['expression_data_CGGA']],
    sep = ','
)

## load the clinical data
clinical_data_CGGA <- fread(
    file = snakemake@input[['clinical_data_CGGA']],
    sep = ','
)

## load the genes of interest
gene_of_interest <- fread(
    file = snakemake@input[['gene_of_interest']]
)

## load the function for creating the heatmap
source(snakemake@params[["heatmap_function"]])
source(snakemake@params[["clustering_metrics"]])
source(snakemake@params[["utils"]])

## generate the directory that will contain all the CGGA validation figures
dir.create(snakemake@output[["cgga_validation_dir"]])


###################################################
## format the data
###################################################


## remove the unknown IDH information
clinical_data_CGGA <- clinical_data_CGGA[
    IDH_mutation_status != 'Unknown',
]

## remove the unknown 1p/19q codeletion information
clinical_data_CGGA <- clinical_data_CGGA[
    `1p19q_codeletion_status` != 'Unknown',
]

## copy the expression data
expression_data_CGGA_filtered <- copy(expression_data_CGGA)

## extract the expression associated with the genes of interest
expression_data_CGGA_filtered <- expression_data_CGGA_filtered[
    Gene_Name %in% unlist(gene_of_interest[, gene_name])
    ,
]

## extract the expression associated with the sample of interest
expression_data_CGGA_filtered <- expression_data_CGGA_filtered[
    ,
    c('Gene_Name', unlist(clinical_data_CGGA[, sample_id])),
    with = F
]


###################################################
## Generate the heatmap of the data
###################################################


## choose the optimal number of cluster with the entropy calculation
#optimal_nb_cluster <- calculate_entropy_varying_clusterNumber(
#    CLUSTER_NUMBER_MAX_INPUT = 20,
#    EXPRESSION_DATA_INPUT = expression_data_CGGA_filtered, # expression input used for the entropy calculation
#    VARIABLE_DATA_INPUT = clinical_data_CGGA, # table containing all the variables data used for the heatmap
#    VARIABLE_INPUT = variable_of_interest, # variables used in the heatmap
#    ID_SAMPLE_VARIABLE_NAMES_INPUT = 'sample_id',
#    PATH_INPUT = snakemake@output[["entropy_curve"]]
#)

## select variables of interest
variable_of_interest <- c(
   'Histology',
   'Grade',
    'IDH_mutation_status',
   '1p19q_codeletion_status'
)

## generate a heatmap for getting the four NT clusters
heatmap <- create_heatmap2(
    expression_data = expression_data_CGGA_filtered,
    clinical = clinical_data_CGGA,
    variables = variable_of_interest,
    genes = unlist(expression_data_CGGA_filtered[, Gene_Name]),
    log_normalization = F,
    col_split = 4,
    clustering_distance_rows = "minkowski",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "spearman",
    clustering_method_columns = "ward.D2"
)


##########
## Extract clusters associated from the heatmap
##########


## load the cluster from the heatmap
cluster_data <- extract_clusters(
    HEATMAP_OBJECT_INPUT = heatmap,
    EXPRESSION_DATA_INPUT = expression_data_CGGA_filtered 
)

## rename each cluster
cluster_data[cluster == 1, cluster := "NT-1 Like"]
cluster_data[cluster == 2, cluster := "NT-2 Like"]
cluster_data[cluster == 3, cluster := "NT-3 Like"]
cluster_data[cluster == 4, cluster := "NT-4 Like"]


##########
## Generate the heatmap and its annotation and save them into file
##########


## merge the cluster data with the clinical data to get the cluster information
clinical_data <- merge(
    clinical_data_CGGA,
    cluster_data,
    by = "sample_id",
    all.x = T,
    sort = F
)

## select variables of interest
variable_of_interest <- c(
    'cluster'
)

## generate a heatmap of the cancer data
heatmap <- create_heatmap2(
    expression_data = expression_data_CGGA_filtered,
    clinical = clinical_data,
    variables = variable_of_interest,
    genes = unlist(expression_data_CGGA_filtered[, Gene_Name]),
    path_file = paste(
        snakemake@output[["cgga_validation_dir"]],
        "/",
        "CGGA_IDHall_heatmap",
        ".svg",
        sep = ""
    ),
    log_normalization = F,
    col_split = 4,
    #row_split = 4,
    clustering_distance_rows = "minkowski",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "spearman",
    clustering_method_columns = "ward.D2"
)

## select variables of interest
variable_of_interest <- c(
    'cluster',
   'Histology',
   'Grade',
    'IDH_mutation_status',
   '1p19q_codeletion_status'
)

## change the annotation for the categories
print(colnames(clinical_data))
clinical_data[Grade == "WHO II", Grade := "G2"]
clinical_data[Grade == "WHO III", Grade := "G3"]
clinical_data[Grade == "WHO IV", Grade := "G4"]

clinical_data[Histology == "A", Histology := "Astrocytoma"]
clinical_data[Histology == "AA", Histology := "Astrocytoma"]
clinical_data[Histology == "AO", Histology := "Oligodendroglioma"]
clinical_data[Histology == "O", Histology := "Oligodendroglioma"]
clinical_data[Histology == "GBM", Histology := "Glioblastoma"]
clinical_data[Histology == "sGBM", Histology := "Glioblastoma"]

clinical_data[clinical_data == "Wildtype"] <- "WT"



### generate a heatmap of the cancer data generating the groups of interest
#heatmap_legend <- create_heatmap_annotation(
#    expression_data = expression_data_CGGA_filtered,
#    clinical = clinical_data,
#    variables = variable_of_interest,
#    genes = unlist(expression_data_CGGA_filtered[, Gene_Name]),
#    path_file = paste(
#        snakemake@output[["cgga_validation_dir"]],
#        "/",
#        "CGGA_IDHall_heatmap_annotation",
#        ".svg",
#        sep = ""
#    ),
#    log_normalization = F,
#    col_split = 4,
#    row_split = 3,
#    clustering_distance_rows = "minkowski",
#    clustering_method_rows = "ward.D2",
#    clustering_distance_columns = "spearman",
#    clustering_method_columns = "ward.D2",
#    horizontal_legend = T
#)














###################################################
## Generate expression profile
###################################################

## transpose the expression data
expression_data_dt <- transpose_datatable(
    expression_data_CGGA_filtered,
    column_name = "Gene_Name",
    "sample_id"

)

## merge with the cluster information
expression_data_dt <- merge(
    cluster_data,
    expression_data_dt,
    by = "sample_id"
)[order(`cluster`),]


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

## rename the list name
names(KEGG_genes_list) <- KEGG_genes_file_names


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
                y = log2(value),
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
                y = log2(value),
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




##########
## generate the expression plot list
##########


## initialize the list that will contain the plots
expression_plot_list <- list()

## initialize the list that will contain the expression data
expression_data_list <- list()

## for each metabolic pathway gene set contained in the KEGG gene list
for (i_gene_set in seq(1, length(KEGG_genes_list))) {
    
    ## retreive the metabolic names
    kegg_name <- names(KEGG_genes_list[i_gene_set])

    ## retreive the genes related to the kegg pathway
    kegg_gene_set <- colnames(expression_data_dt)[colnames(expression_data_dt) %in% unlist(KEGG_genes_list[i_gene_set], use.names = FALSE)]


    ## extract the expression data associated qwith the kegg gene set
    expression_data_kegg <- expression_data_dt[, colnames(expression_data_dt) %in% c("sample_id", "cluster", kegg_gene_set), with = F]
    
    ## melt the data
    expression_data_kegg_melt <- melt.data.table(
        expression_data_kegg,
        id.vars = c("sample_id", "cluster"),
        measure.vars = kegg_gene_set
    )

    ## put the expression subset data into a list for counting 
    expression_data_list[[kegg_name]] <- expression_data_kegg_melt

    ## set colors that will be used for the generation of the plots
    color_palette <- brewer.pal(9, "Set1")

    ## name of the clusters that we want to do the expression profile plot
    cluster_name_vector <- unique(unlist(cluster_data[,cluster]))

    ## create a datatable containing the association of colors for each cluster name
    color_data <- data.table(
        cluster = cluster_name_vector,
        color = color_palette[1:(length(cluster_name_vector))]
    )

    ### set the order of the cluster, putting the Healthy condition at the end
    #cluster_name_vector_vector <- unique(unlist(expression_data_kegg_melt$cluster))
    #cluster_name_vector_vector <- sort(cluster_name_vector_vector[cluster_name_vector_vector != "Healthy"])
    #cluster_name_vector_vector <- c(cluster_name_vector_vector, "Healthy")
    
    #print(cluster_name_vector_vector)
    #break
    ### set the order of the cluster
    #expression_data_kegg_melt$cluster <- factor(
    #    expression_data_kegg_melt$cluster,
    #    levels = cluster_name_vector_vector
    #)

    #print(expression_data_kegg_melt)

    ## do the average expression for each genes
    expression_data_kegg_melt <- expression_data_kegg_melt[, mean(as.numeric(value)), by = c("cluster", "variable"), ][, setnames(.SD, "V1", "value")]

    ## create the average expression profile plot
    plot_all_expression_profile_average <- expression_profile_plot_average(
        expression_data_input = expression_data_kegg_melt,
        cluster_name_vector_input = cluster_name_vector,
        color_table_input = color_data,
        title_name = kegg_name,
        gene_name = FALSE
    )

    ## create the file path that will be used for the saving
    file_path <- paste(
        snakemake@output[['cgga_validation_dir']],
        #"data",
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
        snakemake@output[['cgga_validation_dir']],
        #"data",
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


print("okkkk")
Sys.sleep(1000000000)





break





































###################################################
## Generate percentage for formated data
###################################################


##########
## Formate the data
##########


## add a column to the clinical data related to the IDH and codeletion status
classification_vector <- Map(
    function(x, y) {
        paste(x, y, collapse = "_", sep = "_")
    },
    unlist(clinical_data[, 'IDH_mutation_status', with = F]),
    unlist(clinical_data[, '1p19q_codeletion_status', with = F])
)

## transform the list to vector 
classification_vector <- unlist(classification_vector)

## update the data table that contain the clinical data
clinical_data[, classification := classification_vector]

## change the name of each classification categories
clinical_data[classification == "Mutant_Codel", classification := "IDHmut 1p19q-codel"]
clinical_data[classification == "Mutant_Non-codel", classification := "IDHmut non-1p19q-codel"]
clinical_data[classification == "WT_Codel", classification := "IDHwt"]
clinical_data[classification == "WT_Non-codel", classification := "IDHwt"]


##########
## Call the function on the qualitative variable of interest to get the percentage
##########


## retrieve the variable name
qualitative_variable_vector <- c(
    'Gender',
    'Histology',
    'Grade',
    'classification',
    'IDH_mutation_status',
    '1p19q_codeletion_status'
)

## get the percentage for all the variable
qualitative_percentage_res <- lapply(
    qualitative_variable_vector,
    function(x) {
        get_percentage_qualitative(
            clinical_data,
            x,
            "cluster",
            formated = T
        )
    }
)

## append all the results into a data table
qualitative_percentage_res <- Reduce(
    x = qualitative_percentage_res,
    function(x, y) {rbind(x, y)}
)


##########
## Call the function on the quantitative variable of interest to get the percentage
##########


## retrieve the variable name
quantitative_variable_vector <- c(
    "Age_in_years"
)

## get the percentage for all the variable
quantitative_percentage_res <- lapply(
    quantitative_variable_vector,
    function(x) {
        get_percentage_quantitative(
            clinical_data,
            x,
            "cluster",
            formated = T
        )
    }
)

### append all the results into a data table
#quantitative_statistical_res <- Reduce(
#    x = quantitative_percentage_res,
#    function(x, y) {rbind(x, y)}
#)

## get the data
quantitative_percentage_res <- quantitative_percentage_res[[1]]


##########
## Merge the quantitative and the qualitative statistical results into the same data table
##########


formated_data <- rbind(
    quantitative_percentage_res,
    qualitative_percentage_res
)


######
## create the first line that contain the total of sample for each group
######


## extract the header from the clinical data
header_dt <- as.data.table(table(clinical_data[, cluster]))
header_dt[
    ,
    header := paste(
        V1,
        " (n = ",
        N,
        ")",
        sep = ""
    )
]

## extract the header
header_vector <- c(
    "category",
    unlist(header_dt[, header]),
    "variable"
)
header_dt <- as.data.table(as.list(header_vector))

## rename the colnames of the header
colnames(header_dt) <- colnames(formated_data)

## add into the statistical results
formated_data <- rbindlist(
    list(
        header_dt,
        formated_data
    ),
    use.names = F
)

## write the formated data
fwrite(
    formated_data,
    paste(
        snakemake@output[["cgga_validation_dir"]],
        "/",
        "CGGA_IDHall_formated_data",
        ".csv",
        sep = ""
    ),
    sep = ","
)


######
## Generate the barchart of the classification data
######

classification_data <- get_percentage_qualitative(
    clinical_data,
    'classification',
    "cluster",
    formated = F
)

## set the colors
color_palette <- brewer.pal(9, "Dark2")

## generate the barchart associated with the
classification_piechart <- ggplot(
    data = classification_data,
    aes(
        x = variable,
        y = value,
        fill = category
    )
) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = color_palette) +
    labs(
        x = "CGGA Glioma Cluster",
        y = "Percentage",
    ) +
    theme_classic(base_size = 24)


ggsave(
    paste(
        snakemake@output[["cgga_validation_dir"]],
        "/",
        "CGGA_IDHall_barchart",
        ".svg",
        sep = ""
    ),
    plot = classification_piechart,
    device = "svg",
    width = 12,
    height = 8
)


###################################################
## Epidemiological analysis and statistical analysis
###################################################


##########
## call the function for the statistical analysis for qualitative variable (utils.R)
##########


qualitative_stat_res <- lapply(
    qualitative_variable_vector,
    function(x) {
        statistical_analysis_qualitative(
            clinical_data,
            x,
            'cluster'
        )
    }
)


##########
## call the function for the statistical analysis for quantitative variable (utils.R)
##########


quantitative_stat_res <- lapply(
    quantitative_variable_vector,
    function(x) {
        statistical_analysis_quantitative(
            clinical_data,
            x,
            'cluster'
        )
    }
)


##########
## merge the all the results together
##########

## append the results into a same list
statistical_res_list <- append(
    qualitative_stat_res,
    quantitative_stat_res
)

## append all the results into a data table
statistical_res <- Reduce(
    x = statistical_res_list,
    function(x, y) {rbind(x, y)}
)


##########
## Write the statistical analysis
##########


fwrite(
    statistical_res,
    paste(
        snakemake@output[["cgga_validation_dir"]],
        "/",
        "CGGA_IDHall_epidemio_statistic_res",
        ".csv",
        sep = ""
    ),
    sep = ","
)



###################################################
## Survival analysis
###################################################


## choose the variable that we have to keep for the survival analysis
variable_of_interest <- c(
    'sample_id',
    'OS',
    'vital_status'
)

## merge the cluster data with the clinical data taking the data of the variables of interest
survival_data <- merge(
    x = cluster_data,
    y = clinical_data_CGGA[, variable_of_interest, with = F],
    by = 'sample_id',
    all.x = T
)

## rename each variable names
survival_data <- survival_data[, setnames(
    .SD, 
    c("OS", "vital_status"), 
    c("time", "status"))]


## transform the time to numeric values
survival_data[
    ,
    time := sapply(time, function(x) as.numeric(x))
]

## days_to_death : transform days to months
survival_data[, time := sapply(
    time, 
    function(x) {
        return(as.numeric(x) * 12 / 365)
    }
)]


## transform the vital status to binar
survival_data[, status := ifelse(status == "dead", 1, 0)]

## choose color
color_palette <- brewer.pal(9, "Set1")

## create the survival curve calling the function
create_survival_curve(
    FIT_MODEL_INPUT = survfit(Surv(time, status) ~ cluster, data = survival_data),
    FILE_PATH_INPUT = paste(
        snakemake@output[["cgga_validation_dir"]],
        "/",
        "CGGA_IDHall_survival_curve",
        ".svg",
        sep = ""
    ),
    COLOR_PALETTE_INPUT = color_palette,
    DISPLAY_PVALUE = F,
    HEIGTH_INPUT = 12,
    WIDTH_INPUT = 12
)


##########
## Generate the p values for each cluster comparison for survival
##########

## retreive all the cluster name into a vector
cluster_name_vector <- unique(unlist(cluster_data[, cluster]))

## initialize the data table that will contain the pvalues for each combination
pval_all_dt <- data.table(
    "cluster1" = character(),
    "cluster2" = character(),
    "pval" = numeric()
)

## do the combination of two cluster
cluster_combination <- combn(
    cluster_name_vector,
    m = 2
)

## for each combination, extract the data and generate the survival data
for (i_combination in seq(1, ncol(cluster_combination))) {

    ## retreive the cluster 1
    cluster_name1 <- cluster_combination[1, i_combination]
    cluster_name2 <- cluster_combination[2, i_combination]

    ## copy the subset data of the survival data
    subdata_survival <- copy(survival_data)

    ## extract data related to the combination cluster
    subdata_survival <- subdata_survival[cluster %in% c(cluster_name1, cluster_name2),]

    ## generate the survival data
    sub_fitdata <- survfit(Surv(time, status) ~ cluster, data = subdata_survival)

    ## retreive the pvalues
    pval <- surv_pvalue(sub_fitdata)$pval

    ## put the cluster combination names associated with their pvalue in a data table
    pval_dt <- data.table(
        "cluster1" = cluster_name1,
        "cluster2" = cluster_name2,
        "pval" = pval
    )

    ## rbind the pval_dt with the data table that contain all the pvalue
    pval_all_dt <- rbind(
        pval_all_dt,
        pval_dt
    )
}


## write the data
fwrite(
    pval_all_dt,
    paste(
        snakemake@output[["cgga_validation_dir"]],
        "/",
        "CGGA_IDHall_survival_pvalues",
        ".csv",
        sep = ""
    ),
    sep = ","
)



###################################################
## Split the expression data table for the input of timer2
###################################################


## create the bins
#bins_vector <- c(
#    seq(2, ncol(expression_data_CGGA_filtered), 300),
#    ncol(expression_data_CGGA_filtered)
#)

### create the bins
#for (i_bins in seq(1, (length(bins_vector)-1), 1)) {


#    if (i_bins == (length(bins_vector)-1)) {
#        ## retreive the column number to extract
#        column_nb <- c(
#            1,
#            seq(bins_vector[i_bins], (bins_vector[i_bins + 1]))
#        )

#        ## generate name of the file
#        file_name <- paste(
#            "cgga_expression_",
#            as.character(bins_vector[i_bins]),
#            "_",
#            as.character(bins_vector[i_bins + 1]),
#            ".csv",
#            sep = ""
#        )
#    }

#    else {
#        ## retreive the column number to extract
#        column_nb <- c(
#            1,
#            seq(bins_vector[i_bins], (bins_vector[i_bins + 1]-1))
#        )

#        ## generate name of the file
#        file_name <- paste(
#            "cgga_expression_",
#            as.character(bins_vector[i_bins]),
#            "_",
#            as.character(bins_vector[i_bins + 1]-1),
#            ".csv",
#            sep = ""
#        )
#    }

#    ## generate the path of the file
#    file_path <- paste(
#        snakemake@output[["cgga_validation_dir"]],
#        "/",
#        file_name,
#        sep = ""
#    )

#    ## retrive the expression data associated with the column number
#    subset_data <- copy(expression_data_CGGA_filtered[,column_nb, with = F])

#    ## transform the data table to data frame for saving for TIMER2
#    setDF(subset_data, rownames = unlist(subset_data[,Gene_Name]))

#    ## drop column from the data frame 
#    drop_col <- c("Gene_Name")

#    ## remove the column that contain the Gene_Name
#    subset_data <- subset_data[, !(colnames(subset_data) %in% drop_col)]

#    write.csv(
#        subset_data,
#        file_path,
#        sep = ",",
#        quote = FALSE
#    )

#}
#Sys.sleep(10000)

###################################################
## Inferring cell type 
###################################################

## load the files contained in the directory related to the cgga immune estimation
immune_file_path <- list.files(
    snakemake@input[["cgga_immune_estimation_dir"]],
    full.names = T
)

## load each file
immune_data_list <- lapply(
    immune_file_path,
    function(x) fread(
        x,
        sep = ","
    )
)

## merge all the data into a same data table
immune_data_dt <- Reduce(
    x = immune_data_list,
    function(x, y) {
        merged_data <- merge(
            x,
            y,
            by = "cell_type",
            all.x = T,
            sort = F
        )
    }
)

## transpose the immune data table
immune_data_dt <- transpose(
    immune_data_dt,
    keep.names = "sample_id",
    make.names = "cell_type"

)

## merge the immune data with the cluster data
immune_data_dt <- merge(
    cluster_data,
    immune_data_dt,
    by = "sample_id",
    all.x = T,
    sort = F
)

## color to use in the plot
all_colors <- c(
    brewer.pal(9, "Set1"),
    brewer.pal(8, "Dark2"),
    brewer.pal(8, "Set2"),
    brewer.pal(12, "Paired"),
    brewer.pal(9, "Set1")
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
    colname_of_interest <- grep(x = colnames(immune_data_dt), pattern = tools_regex, value = T)

    ## extract the data associated with the column name of interest
    subset_data <- immune_data_dt[,c('sample_id', 'cluster', colname_of_interest), with = F]
    subset_data <- immune_data_dt[,c('cluster', colname_of_interest), with = F]

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
        '.svg',
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
        device = 'svg',
        width = 18,
        height = 8
    )

    return(immune_plot)
}


## set the data of which tools we will used
tools <- c(
    'CIBERSORT',
    'CIBERSORT-ABS',
    'QUANTISEQ'
#    'EPIC',
#    'TIMER',
#    'XCELL'
)

## create the immune plots for every selected tools
immune_plot_list <- lapply(
    tools,
    function(x) create_immune_stacked_barplot(
        x,
        directory_path = snakemake@output[["cgga_validation_dir"]])
)


##########
## Do the statistical analysis between each clusters
##########

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


## extract all the data related to the tools of interest into a list
immune_data_list <- lapply(
    tools,
    function(x) extract_immune_data(
        tool_name = x,
        immune_data = immune_data_dt
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
            snakemake@output[["cgga_validation_dir"]],
            "/",
            names(statistical_results_list[i_stat_results]),
            "_pval",
            ".csv",
            sep = ""
    )
    
    print(statistical_results_list[[i_stat_results]])
    ## write the data table that contain the stats into a file
    fwrite(
        x = (statistical_results_list[[i_stat_results]]),
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
    immune_data_dt

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
            snakemake@output[["cgga_validation_dir"]],
            "/",
            "pval_cibersortabs_quantification_per_cluster",
            ".csv",
            sep = ""
    ),
    sep = ","
)
