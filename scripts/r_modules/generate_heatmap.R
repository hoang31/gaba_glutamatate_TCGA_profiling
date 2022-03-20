
###################################################
###################################################


## Create heatmap from the expressiom data


###################################################
###################################################


## load packages
library("data.table")
library("RColorBrewer")
library("ComplexHeatmap")
suppressMessages(library("tidyverse"))

#library("NbClust")

###################################################
###################################################


##### load data

## load emsebl and symbols id associated with the GABA, GLUTA and CALCIUM pathways
gene_id_of_interest <- fread(
    snakemake@input[["gene_id"]],
    sep = ",",
    header = T
)

## load expression data
d_expression <- fread(
    file = snakemake@input[["expression_data"]],
    sep = ",",
    header = T
)

## load clinical data
d_clinical_TCGA_normal <- fread(
    input = snakemake@input[["TCGA_clinical_normal"]],
    sep = ",",
    header = T
)

d_clinical_TCGA_cancer <- fread(
    input = snakemake@input[["TCGA_clinical_cancer"]],
    sep = ",",
    header = T
)

d_clinical_TCGA_external <- fread(
    input = snakemake@input[["TCGA_clinical_EXTERNAL"]],
    sep = ",",
    header = T
)

## load external scripts
source(snakemake@params[["heatmap_function"]])
source(snakemake@params[["clustering_metrics"]])
source(snakemake@params[["utils_R"]])


###################################################
###################################################

##### EXTRACT THE EXPRESSION DATA OF THE GENES OF INTEREST

## tranform the ensembl id to gene symbol and extract the gene of interest from the gene_id data
d_expression[, rn := ensemblID_to_geneSymbole(rn, gene_id_of_interest)]
d_expression <- d_expression[!grep(rn, pattern = "ENS"),][, setnames(.SD, "rn", "Gene_Name")]

##### PARSE THE CLINICAL DATA

## rbind the TCGA data associated with cancer data and normal data
d_clinical_TCGA <- rbind(
    d_clinical_TCGA_normal,
    d_clinical_TCGA_cancer
)

## extract the id for the count files
d_clinical_TCGA <- unique(
    d_clinical_TCGA[
        grep(
            x = d_clinical_TCGA[, File_Name],
            pattern = "counts"
        ),
    ]
)

## create a column called "type" associated with the tumor type, normal or cancer
d_clinical_TCGA[
    ,
    "type" := sapply(
        d_clinical_TCGA[, Sample_Type],
        function(x) {
            ifelse(
                x == "Solid Tissue Normal",
                "normal",
                "cancer"
            )
        }
    )
]

## remove duplicate column (Project ID)
d_clinical_TCGA <- d_clinical_TCGA[
    ,
    unique(names(d_clinical_TCGA)),
    with = FALSE
]

## merge clinical data with other clinical data
d_clinical_TCGA <- merge(
    d_clinical_TCGA,
    d_clinical_TCGA_external,
    by.x = "Case_ID",
    by.y = "Case",
    all.x = T,
    sort = F
)

## create a column corresponding to the database name
d_clinical_TCGA[, "database" := "TCGA"]

## NORMAL LIKE glioma previously identified
NL_samples <- c(
    "TCGA-DU-A76K",
    "TCGA-P5-A5EY",
    "TCGA-DB-A75P",
    "TCGA-HT-8558",
    "TCGA-CS-6669",
    "TCGA-DU-7292",
    "TCGA-QH-A6XC",
    "TCGA-HT-8015",
    "TCGA-HT-8564",
    "TCGA-FG-8181",
    "TCGA-HT-8107",
    "TCGA-HT-8019",
    "TCGA-DU-8162",
    "TCGA-FG-7643"
)

d_clinical_TCGA[Case_ID %in% NL_samples, "cancer type" := "NL" ]
d_clinical_TCGA[!(Case_ID %in% NL_samples), "cancer type" := "OT" ]

## create a variable called "classification" describing the glioma classification : IDH-mut codel, IDH-mut noncodel and IDHwt
d_clinical_TCGA[`IDH status` == "WT", "classification" := "IDHwt"]
d_clinical_TCGA[(`IDH status` == "Mutant") & (`1p/19q codeletion` == "codel"), "classification" := "IDHmut_CODEL"]
d_clinical_TCGA[(`IDH status` == "Mutant") & (`1p/19q codeletion` == "non-codel"), "classification" := "IDHmut_nonCODEL"]




####################################################
####################################################


##### SEPARATION OF IDH-WT, IDH-MUT with 1p/19q codel and IDH-mut without the 1p/19q codel


## extract samples associated with each class
d_clinical_IDHwt <- copy(d_clinical_TCGA[(classification == "IDHwt"), ])

## IDH mutated
d_clinical_IDHmut <- copy(d_clinical_TCGA[(classification == "IDHmut_CODEL") | (classification == "IDHmut_nonCODEL"), ])
d_clinical_IDHmut_CODEL <- copy(d_clinical_TCGA[(classification == "IDHmut_CODEL"), ])
d_clinical_IDHmut_nonCODEL <- copy(d_clinical_TCGA[(classification == "IDHmut_nonCODEL"), ])

## extract the sample assoiated with the glioma samples
d_clinical_cancer <- copy(d_clinical_TCGA[(type == "cancer") & classification %in% c("IDHwt", "IDHmut_CODEL", "IDHmut_nonCODEL"), ])








## remove the glioblastoma
#d_clinical_cancer <- d_clinical_cancer[Histology != 'glioblastoma',]









## extract expression data associated with the IDH-wt and IDH-mut cases
d_expression_IDHwt <- copy(d_expression[
    ,
    append(
        "Gene_Name",
        unlist(d_clinical_IDHwt[, "Case_ID"])
    ),
    with = F
])

d_expression_IDHmut <- copy(d_expression[
    ,
    append(
        "Gene_Name",
        unlist(d_clinical_IDHmut[, "Case_ID"])
    ),
    with = F
])

d_expression_IDHmut_CODEL <- copy(d_expression[
    ,
    append(
        "Gene_Name",
        unlist(d_clinical_IDHmut_CODEL[, "Case_ID"])
    ),
    with = F
])

d_expression_IDHmut_nonCODEL <- copy(d_expression[
    ,
    append(
        "Gene_Name",
        unlist(d_clinical_IDHmut_nonCODEL[, "Case_ID"])
    ),
    with = F
])


d_expression_cancer <- copy(d_expression[
    ,
    append(
        "Gene_Name",
        unlist(d_clinical_cancer[, "Case_ID"])
    ),
    with = F
])

cat("\n###################################################\n")
cat("#  ----- TCGA LGG GBB Summary", "-----", "\n")
cat("#   Total cancer samples :", ncol(d_expression_cancer) - 1, "\n")
cat("#   IDH-wt :", ncol(d_expression_IDHwt) - 1, "\n")
cat("#   IDH-mut :", ncol(d_expression_IDHmut) - 1, "\n")
cat("#     -> IDH-mut with 1p/19q codeletion :", ncol(d_expression_IDHmut_CODEL) - 1, "\n")
cat("#     -> IDH-mut without 1p/19q codeletion :", ncol(d_expression_IDHmut_nonCODEL) - 1, "\n")
cat("###################################################\n\n")


###################################################
###################################################


##### GENERATE THE HEATMAPS 


###################################################


## initialize the variables
# var <- c(
#     "type",
#     "IDH status",
#     "Histology",
#     "Grade",
#     "1p/19q codeletion",
#     "Chr 7 gain/Chr 10 loss",
#     "TERT promoter status",
#     "ATRX status",
#     "BRAF V600E status",
#     "vital_status"
# )


##### IDH ALL (mutated + wild type)

## initialize the variables
var <- c(
    "IDH status",
    "1p/19q codeletion"
)

optimal_nb_cluster_IDHall <- calculate_entropy_varying_clusterNumber(
    CLUSTER_NUMBER_MAX_INPUT = 20,
    EXPRESSION_DATA_INPUT = d_expression_cancer, # expression input used for the entropy calculation
    VARIABLE_DATA_INPUT = d_clinical_cancer, # table containing all the variables data used for the heatmap
    VARIABLE_INPUT = var, # variables used in the heatmap
    ID_SAMPLE_VARIABLE_NAMES_INPUT = 'Case_ID',
    PATH_INPUT = snakemake@output[["optimal_nb_clusters_IDHall"]]
)






















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

## extract all unique genes
unique_kegg_genes <- unique(unlist(KEGG_genes_list, recursive = T))

## filter the kegg genes with the expressed genes
unique_kegg_genes <- unique_kegg_genes[unique_kegg_genes %in% unlist(d_expression_cancer[, Gene_Name])]

## create a data table that contains these geness
unique_kegg_genes_dt <- data.table(
    gene_name = unique_kegg_genes
)

## for each kegg genes, put the information
for (i_kegg_list in seq(1, length(KEGG_genes_file_names))) {

    ## extract the kegg name and the genes associated with this latter
    kegg_gene_vector <- KEGG_genes_list[[i_kegg_list]]
    kegg_name <- names(KEGG_genes_list[i_kegg_list])
    
    ## rename the genes
    unique_kegg_genes_dt[gene_name %in% kegg_gene_vector, kegg_pathway := kegg_name]
    print(kegg_gene_vector)
}






## function for creating the annotation for the heatmap
create_row_annotation <- function(
    variable,
    clinicalData,
    column_annotation = TRUE,
    color_vector
) {
    ## create the colors
    #color_vector <- c(
    #    brewer.pal(9, "Set1"),
    #    brewer.pal(8, "Dark2"),
    #    brewer.pal(8, "Set2"),
    #    brewer.pal(12, "Paired")
    #)

    ## extract data table with the specific variables
    clinicalData <- clinicalData[, variable, with = F]

    ## initialise the color list
    colors_for_variables <- list()

    ## for each data associated with each variables, do a counting of occurence and extract colors
    for (i in (1:length(variable)))
    {

        ## extract the number of catagories in the variable i-th
        categories <- names(table(clinicalData[, variable[i], with = F]))

        ## take the i-th first colors from the "color_vector" variable
        colors_for_categories <- color_vector[1:length(categories)]

        ## remove the colors used from the initial color vector
        color_vector <- color_vector[!(color_vector %in% colors_for_categories)]
 

        ## name each color by the category
        names(colors_for_categories) <- categories

        ## append into the list the colors categories
        colors_for_variables[[i]] <- colors_for_categories
    }

    ## give variable names for each vector in the list
    names(colors_for_variables) <- variable


 
    ## create the annotations
    annotation <- rowAnnotation(
        df = clinicalData,
        col = colors_for_variables
    )
    return(annotation)


}


##### function for creating a heatmap
create_heatmap2 <- function(expression_data,
                            clinical,
                            variables,
                            genes,
                            pseudo_count,
                            path_file,
                            col_split,
                            row_split,
                            return_matrix,
                            log_normalization = T,
                            clustering_distance_rows = "euclidean",
                            clustering_method_rows = "complete",
                            clustering_distance_columns = "euclidean",
                            clustering_method_columns = "complete") {

    ## generate the heatmap file
    generation_heatmap <- FALSE

    ## return the count matrix ?
    if (missing(return_matrix)) {
        return_matrix <- FALSE
    }

    if (!(missing(path_file))) {
        generation_heatmap <- TRUE
    }

    ## boolean for splitting or not the rows and columns with the clustering
    split_col <- FALSE
    split_row <- FALSE

    ## if the argument col_split or row_split are present
    if (!(missing(col_split))) {
        split_col <- TRUE
    }

    if (!(missing(row_split))) {
        split_row <- TRUE
    }

    ## if there is a column spliting
    if (col_split == 0) {
        split_col <- FALSE
    }

    if (split_row == 0) {
        split_row <- FALSE
    }

    ## extract the expression related to the specific genes
    expression_filtered <- expression_data[Gene_Name %in% genes, ]

    ## transform the data table to matrix removing the "Gene_Name" column
    expression_filtered_matrix <- as.matrix(
        x = expression_filtered[, !c("Gene_Name")],
        rownames = expression_filtered[, Gene_Name]
    )

    ## if we return the filtered matrix
    if (return_matrix == TRUE) {
        return(expression_filtered_matrix)
    }

    if (log_normalization == TRUE) {
        # print("log normalization...")
        ## transform the matrix with a log + pseudocount
        expression_filtered_matrix_log <- apply(expression_filtered_matrix, c(1, 2), function(x) log2(x + pseudo_count))
    }

    if (log_normalization == FALSE) {
        # print("no log normalization...")
        ## transform the matrix with a log + pseudocount
        expression_filtered_matrix_log <- copy(expression_filtered_matrix)
    }

    ## extract and make the annotations matrix
    annotation <- create_annotation(variables, clinical)


    gene_information <- merge(
        expression_filtered,
        unique_kegg_genes_dt,
        by.x = "Gene_Name",
        by.y = "gene_name",
        all.x = T,
        sort = F
    )

    ## create the row annotation
    row_annotations <- create_row_annotation(
        variable = c(
            'kegg_pathway'
        ),
        clinicalData = gene_information,
        color_vector = brewer.pal(8, "Set2")
    )


    ## write the file
    if (generation_heatmap) { svg(path_file, width = 10, height = 8) } 

    ## color
    ht_colors <- colorRamp2(
        c(6, 9, 14),
        c("blue3", "white", "red3")
    )


    column_title <- "SAMPLES"
    row_title <- "GENES"
    column_title_gp <- gpar(fontsize = 22, fontface = "bold")
    row_title_gp <- gpar(fontsize = 22, fontface = "bold")

    ## make the heatmap with the row and column splitting
    if ((split_col == T) & (split_row == T)) {
        # print("split rows and columns")
        ht <- (Heatmap(expression_filtered_matrix_log,
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,

            column_gap = unit(5, "mm"),
            row_gap = unit(5, "mm"),

            border = TRUE,
            show_column_names = FALSE,
            show_row_names = FALSE,
            column_split = col_split,
            row_split = row_split,

            top_annotation = annotation,
            left_annotation = row_annotations,

            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),

            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,

            row_dend_reorder = FALSE,
            column_dend_reorder = TRUE
        )
        )
    }

    ## make the heatmap with the column splitting
    if ((split_col == T) & (split_row == F)) {
        # print("split columns")
        ht <- (Heatmap((expression_filtered_matrix_log),
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,
            border = TRUE,
            show_column_names = FALSE,
            show_row_names = FALSE,
            column_split = col_split,
            top_annotation = annotation,
            left_annotation = row_annotations,
            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),
            column_gap = unit(5, "mm"),
            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,

            row_dend_reorder = FALSE,
            column_dend_reorder = TRUE
        )
        )
    }

    ## make the heatmap with the row splitting
    if ((split_col == F) & (split_row == T)) {
        # print("split rows")
        ht <- (Heatmap((expression_filtered_matrix_log),
            border = TRUE,
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,
            show_column_names = FALSE,
            show_row_names = FALSE,
            row_split = row_split,
            row_gap = unit(5, "mm"),
            top_annotation = annotation,
            left_annotation = row_annotations,
            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),
            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,
            row_dend_reorder = FALSE,
            column_dend_reorder = TRUE
        )
        )
    }

    ## make the heatmap without the splitting
    if ((split_col == F) & (split_row == F)) {
        # print("no split rows and columns")
        ht <- (Heatmap((expression_filtered_matrix_log),
            border = TRUE,
            col = ht_colors,
            column_title = column_title,
            row_title = row_title,
            column_title_gp = column_title_gp,
            row_title_gp = row_title_gp,
            show_column_names = FALSE,
            show_row_names = FALSE,
            ## legends for the heatmap
            top_annotation = annotation,
            left_annotation = row_annotations,
            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                at = c(0, 5, 10, 15),
                grid_width = unit(0.75, "cm"),
                legend_height = unit(6, "cm"),
                border = "black"
            ),
            clustering_distance_rows = clustering_distance_rows,
            clustering_method_rows = clustering_method_rows,
            clustering_distance_columns = clustering_distance_columns,
            clustering_method_columns = clustering_method_columns,
            row_dend_reorder = F,
            column_dend_reorder = F
        )
        )
    }


    if (generation_heatmap) {
        print(ht)
        dev.off()
    }

    return(ht)

    print("###############################################################")
}





























































































































## generate a heatmap of the cancer data without generate the groups
heatmap_IDHall <- create_heatmap2(
    expression_data = d_expression_cancer,
    clinical = d_clinical_cancer,
    variables = var,
    genes = unlist(d_expression_cancer[, Gene_Name]),
    pseudo_count = 0.1,
    path_file = "data/figures/TCGA_heatmap_without_groups.svg",
    log_normalization = F,
    col_split = 0,
    clustering_distance_rows = "minkowski",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "spearman",
    clustering_method_columns = "ward.D2"
)




## generate a heatmap of the cancer data generating the groups of interest
heatmap_IDHall <- create_heatmap2(
    expression_data = d_expression_cancer,
    clinical = d_clinical_cancer,
    variables = var,
    genes = unlist(d_expression_cancer[, Gene_Name]),
    pseudo_count = 0.1,
    path_file = snakemake@output[["heatmap_IDHall"]],
    log_normalization = F,
    col_split = 4,
    row_split = 3,
    clustering_distance_rows = "minkowski",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "spearman",
    clustering_method_columns = "ward.D2"
)




## heatmap with the normal like informations
var <- c(
    "IDH status",
    "1p/19q codeletion",
    "cancer type"
)

heatmap_NL <- create_heatmap2(
    expression_data = d_expression_cancer,
    clinical = d_clinical_cancer,
    variables = var,
    genes = unlist(d_expression_cancer[, Gene_Name]),
    pseudo_count = 0.1,
    path_file = "data/figures/TCGA_heatmap_with_NLsamples.svg",
    log_normalization = F,
    col_split = 4,
    clustering_distance_rows = "minkowski",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "spearman",
    clustering_method_columns = "ward.D2"
)


# ## if we generate four clusters
# if (optimal_nb_cluster_IDHall == 4) {


#     ###### reorder the clusters

#     ## get the colomn clusters
#     column_clusters <- column_order(heatmap_IDHall)
#     samples_ids <- data.table(
#         sample_id = colnames(d_expression_cancer)[colnames(d_expression_cancer) != "Gene_Name"]
#     )

#     ## get the cluster names
#     cluster_names <- c(
#         "MIXED",
#         "IDH_WT",
#         "IDH_MUT_codel",
#         "IDH_MUT_noncodel"
#     )

#     ## rename the cluster by the cluster name
#     for (i_column_cluster in seq(1, length(column_clusters), 1)) {
#         pos <- column_clusters[[i_column_cluster]]
#         samples_ids[pos, cluster:= cluster_names[i_column_cluster]]
#     }

#     ## transform as a vector
#     cluster_vector <- factor(
#         unlist(samples_ids[, "cluster"]),
#         levels = c(
#             "IDH_MUT_codel",
#             "IDH_MUT_noncodel",
#             "IDH_WT",
#             "MIXED"
#         )
#     )

#     table(cluster_vector)

#     ## generate the heatmap
#     heatmap_IDHall2 <- create_heatmap3(
#         expression_data = d_expression_cancer,
#         clinical = d_clinical_cancer,
#         variables = var,
#         genes = unlist(d_expression_cancer[, Gene_Name]),
#         pseudo_count = 0.1,
#         path_file = "data/figures/heatmap_ordered.svg",
#         log_normalization = F,
#         clustering_distance_rows = "minkowski",
#         clustering_method_rows = "ward.D2",
#         clustering_distance_columns = "spearman",
#         clustering_method_columns = "ward.D2",
#         col_split = 0,
#         cluster_name = cluster_vector
#     )

# }




###################################################


##### IDH WT

# ind <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "dunn", "sdindex", "sdbw")

## search the optimal number of clusters : compute for all the index with nbclust
# optimal_nb_cluster_IDHwt <- calculate_nbclust_for_indexes(
#     EXPRESSION_DATA_INPUT = d_expression_IDHwt,
#     SUBTITLE_INPUT = "IDH-WT glioma - samples",
#     INDEX_VECTOR_INPUT = ind,
#     WHERE_INPUT = "column",
#     DISTANCE_INPUT = "spearman",
#     PATH_INPUT = snakemake@output[["optimal_nb_clusters_IDHwt"]]
# )

## initialize the variables
# var <- c(
    # "IDH status",
    # "Histology",
    # "Grade",
    # "vital_status"
# )

## generate a heatmap of the IDH-wt data
# heatmap_IDHwt <- create_heatmap2(
#     expression_data = d_expression_IDHwt,
#     clinical = d_clinical_IDHwt,
#     variables = var,
#     genes = unlist(d_expression_IDHwt[, Gene_Name]),
#     pseudo_count = 0.1,
#     path_file = snakemake@output[["heatmap_IDHwt"]],
#     log_normalization = F,
#     col_split = 2,
#     clustering_distance_rows = "minkowski",
#     clustering_method_rows = "ward.D2",
#     clustering_distance_columns = "spearman",
#     clustering_method_columns = "ward.D2"
# )


###################################################


##### IDH MUTATED

# ind <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "dunn", "sdindex", "sdbw")

# search the optimal number of clusters : compute for all the index with nbclust
# optimal_nb_cluster_IDHmut <- calculate_nbclust_for_indexes(
#     EXPRESSION_DATA_INPUT = d_expression_IDHmut,
#     SUBTITLE_INPUT = "IDH-MUT glioma - samples",
#     INDEX_VECTOR_INPUT = ind,
#     WHERE_INPUT = "column",
#     DISTANCE_INPUT = "spearman",
#     PATH_INPUT = snakemake@output[["optimal_nb_clusters_IDHmut"]]
# )

# var <- c(
#     "IDH status",
#     "1p/19q codeletion",
#     "Histology",
#     "Grade",
#     "vital_status"
# )

## generate a heatmap of the IDH-mut data
# heatmap_IDHmut <- create_heatmap2(
#     expression_data = d_expression_IDHmut,
#     clinical = d_clinical_IDHmut,
#     variables = var,
#     genes = unlist(d_expression_IDHmut[, Gene_Name]),
#     pseudo_count = 0.1,
#     path_file = snakemake@output[["heatmap_IDHmut"]],
#     log_normalization = F,
#     col_split = 2,
#     clustering_distance_rows = "minkowski",
#     clustering_method_rows = "ward.D2",
#     clustering_distance_columns = "spearman",
#     clustering_method_columns = "ward.D2"
# )


###################################################
###################################################


##### EXTRACT THE CLUSTERING GROUPS
## function for extracting the clustering groups
extract_clusters <- function(
    HEATMAP_OBJECT_INPUT, # object generated with the the function "heatmap2"
    EXPRESSION_DATA_INPUT # expression data used for the generation of the heatmap
) {

    ## extract the clustering groups of the heatmap
    clustering_groups <- column_order(HEATMAP_OBJECT_INPUT)

    ## extract the samples from the expression data
    samples_groups <- as.data.table(
        list(colnames(EXPRESSION_DATA_INPUT)[colnames(EXPRESSION_DATA_INPUT) != "Gene_Name"])
    )[
        ,
        setnames(
            .SD,
            "V1",
            "sample_id"
        )
    ]

    ## for each clustering groups, add the group number to the "samples_groups" data frames
    for (i in seq(1, length(clustering_groups), 1)) {
        # clustering_groups[clustering_groups[[i]], cluster := i]

        set(
            samples_groups,
            i = clustering_groups[[i]],
            j = "cluster",
            value = toString(i)
        )
    }

    return(samples_groups)

}

## extract cluster groups of each generated heatmap
cluster_IDHall <- extract_clusters(
    heatmap_IDHall,
    d_expression_cancer
)


###################################################################
###################################################################


## for renaming clustering depending of the number of glioma type that compose the group

#clustering_data <- merge(
#    cluster_IDHall,
#    d_clinical_TCGA[,c('submitter_id', 'classification')],
#    by.x = "sample_id",
#    by.y = "submitter_id",
#)

### counting for each glioma type of each cluster 
#classification_counting <- as.data.table(as.data.frame.matrix(table(clustering_data[, c('cluster', 'classification')])), keep.rownames = "cluster_name")

### do the sum by cluster
#classification_counting[, total := sum(.SD), by = 'cluster_name']


#classification_counting_percentage <- classification_counting[
#    ,
#    `:=` (
#        'IDHmut_CODEL_percentage' = IDHmut_CODEL/total*100,
#        'IDHmut_nonCODEL_percentage' = IDHmut_nonCODEL/total*100,
#        'IDHwt_percentage' = IDHwt/total*100
#    ),
#    by = 'cluster_name'
#]

#classification_counting_percentage <- classification_counting_percentage[, c('cluster_name', 'IDHmut_CODEL_percentage', 'IDHmut_nonCODEL_percentage', 'IDHwt_percentage'), with = F]

#for (i in seq(1, nrow(classification_counting_percentage), 1)) {
    
#    ## extract the percentage for each glioma type
#    percentage_vector <- unlist(classification_counting_percentage[
#        i,
#        c(
#            'IDHmut_CODEL_percentage',
#            'IDHmut_nonCODEL_percentage',
#            'IDHwt_percentage'
#        ),
#        with = F
#    ])

#    ## extract the percentage that are the highest 
#    highest_percentage <- percentage_vector[order(percentage_vector)[length(order(percentage_vector))]]

#    ## if the percentage is higher than 50%, rename the cluster name by the glioma type names
#    if (as.numeric(highest_percentage) > 50) {

#        ## extract the name of the glioma type that is associated with the highest percentage
#        new_cluster_name <- unlist(strsplit(names(highest_percentage), split = '[_]')[1])
#        new_cluster_name <- new_cluster_name[-(length(new_cluster_name))]
#        new_cluster_name <- paste(new_cluster_name, sep = '', collapse = '_')
#    }

#    ## if all the glioma type percentage are above than 50%, it seems to have all type of glioma so we will name this cluster as a MIXED cluster
#    else {
#        new_cluster_name <- 'MIXED'
#    }

#    ## put the new names associated with the same old name into a same data table
#    classification_counting_percentage[i, cluster_name2 := new_cluster_name]
#}

### put in a independant data table the old and the new cluster name
#cluster_names <- copy(classification_counting_percentage[, c('cluster_name', 'cluster_name2'), with = F])

### extract into vectors the new and the old cluster names
#old_cluster_names <- unique(unlist(cluster_names[, cluster_name])) 
#new_cluster_names <- unique(unlist(cluster_names[, cluster_name2])) 

### if the length of the new cluster name vector is not the same of the length of the old cluster name vector, create a error
#if (length(old_cluster_names) != length(new_cluster_names)) {
#    stop('Error : after renaming the cluster names, some clusters are composed of the same mojority type.')
#}

### function for extrating the new cluster name from the old cluster name
#search_new_cluster_name <- function(old_cluster_name_input, cluster_names_datatable_input) {
#    return(unlist(cluster_names_datatable_input[cluster_name == old_cluster_name_input, cluster_name2]))
#}

### rename the cluster name by the new generated names
#cluster_IDHall[
#    ,
#    cluster := sapply(
#        cluster, 
#        function(x) search_new_cluster_name(x, cluster_names)
#    )
#]




cluster_IDHall[, cluster := paste("NT-", cluster, sep = "")]


## rename the cluster names mannually
# cluster_IDHall[cluster == 1, cluster := "MIXED"]
# cluster_IDHall[cluster == 2, cluster := "IDH_WT"]
# cluster_IDHall[cluster == 3, cluster := "IDH_MUT_codel"]
# cluster_IDHall[cluster == 4, cluster := "IDH_MUT_noncodel"]


###################################################################
###################################################################


## for each cluster, we will create subgroups associated with the histological
cluster_IDHall_histology <- merge(
    x = cluster_IDHall,
    y = d_clinical_TCGA[, c("submitter_id", "primary_diagnosis","Grade", "classification", "1p/19q codeletion"), with = F],
    by.x = "sample_id",
    by.y = "submitter_id",
    all.x = T
)


###################################################
###################################################


#### create subgroups of  from the histological and biomolecular data

### rename the histological entities in the 'primary_diagnosis' column
#cluster_IDHall_histology_biomolecular <- copy(cluster_IDHall_histology)
## cluster_IDHall_histology_biomolecular[, cluster2 := primary_diagnosis]
#cluster_IDHall_histology_biomolecular[primary_diagnosis == 'Astrocytoma, NOS', primary_diagnosis := 'astrocytoma']
#cluster_IDHall_histology_biomolecular[primary_diagnosis == 'Astrocytoma, anaplastic', primary_diagnosis := 'astrocytoma']
#cluster_IDHall_histology_biomolecular[primary_diagnosis == 'Mixed glioma', primary_diagnosis := 'astrocytoma']
#cluster_IDHall_histology_biomolecular[primary_diagnosis == 'Glioblastoma', primary_diagnosis := 'glioblastoma']

### oligodendroglioma
#cluster_IDHall_histology_biomolecular[classification == 'IDHmut_CODEL', cluster2 := 'IDHmut_1p19q_oligodendroglioma']

### glioblastoma
#cluster_IDHall_histology_biomolecular[(classification == 'IDHwt') & (primary_diagnosis == 'glioblastoma'), cluster2 := 'IDHwt_glioblastoma']
#cluster_IDHall_histology_biomolecular[(classification == 'IDHmut_nonCODEL') & (primary_diagnosis == 'glioblastoma'), cluster2 := 'IDHmut_glioblastoma']

### astrocytoma
#cluster_IDHall_histology_biomolecular[(classification == 'IDHmut_nonCODEL') & (primary_diagnosis == 'astrocytoma'), cluster2 := 'IDHmut_astrocytoma']
#cluster_IDHall_histology_biomolecular[(classification == 'IDHmut_nonCODEL') & (primary_diagnosis == 'astrocytoma') & (Grade == 'G2'), cluster2 := 'IDHmut_astrocytoma']
#cluster_IDHall_histology_biomolecular[(classification == 'IDHmut_nonCODEL') & (primary_diagnosis == 'astrocytoma') & (Grade == 'G3'), cluster2 := 'IDHmut_astrocytoma']

#cluster_IDHall_histology_biomolecular[(classification == 'IDHwt') & (primary_diagnosis == 'astrocytoma'), cluster2 := 'IDHwt_astrocytoma']
#cluster_IDHall_histology_biomolecular[(classification == 'IDHwt') & (primary_diagnosis == 'astrocytoma') & (Grade == 'G2'), cluster2 := 'IDHwt_astrocytoma']
#cluster_IDHall_histology_biomolecular[(classification == 'IDHwt') & (primary_diagnosis == 'astrocytoma') & (Grade == 'G3'), cluster2 := 'IDHwt_astrocytoma']

### put unknown for all the na values
#cluster_IDHall_histology_biomolecular[is.na(cluster_IDHall_histology_biomolecular),] <- 'unknown'

### tranform the glioma with mutation in IDH or no mutation and without the codeletion 1p/19q to astrocytoma IDHmut or astrocytoma IDHwt
#cluster_IDHall_histology_biomolecular[(classification == 'IDHmut_nonCODEL') & (cluster2 == 'unknown'), cluster2 := 'IDHmut_astrocytoma']
#cluster_IDHall_histology_biomolecular[(classification == 'IDHwt') & (cluster2 == 'unknown'), cluster2 := 'IDHwt_astrocytoma']

### remove the variable that will be not usefull for the next analysis
#cluster_IDHall_histology_biomolecular <- cluster_IDHall_histology_biomolecular[, c('sample_id', 'cluster','cluster2')]

### add the cluster names into the histological and biomolecular entity names
#cluster_IDHall_histology_biomolecular[
#    ,
#    cluster3 := paste(
#        cluster,
#        '_cluster_',
#        cluster2,
#        sep = ''
#    )
#]

### extract the column named 'cluster3'
#cluster_IDHall_histology_biomolecular <- cluster_IDHall_histology_biomolecular[, c('sample_id', 'cluster3')]

### rename the column
#cluster_IDHall_histology_biomolecular <- cluster_IDHall_histology_biomolecular[
#    ,
#    setnames(
#        .SD,
#        c('cluster3'),
#        c('cluster')
#    )
#]

###################################################
###################################################


## write the data

## write the subcluster data
#fwrite(
#    x = cluster_IDHall_histology_biomolecular,
#    file = snakemake@output[["cluster_IDHall_subclusters"]],
#    sep = ","
#)

# ## extract the cluster position from the heatmap
# cluster_IDHwt <- extract_clusters(
#     heatmap_IDHwt,
#     d_expression_IDHwt
# )

# ## extract the cluster position from the heatmap
# cluster_IDHmut <- extract_clusters(
#     heatmap_IDHmut,
#     d_expression_IDHmut
# )

## put the cluster rresults and the output path in a list
clusters_results_list <- list(
    list(cluster_IDHall, snakemake@output[["cluster_IDHall"]])
    # list(cluster_IDHwt, snakemake@output[["cluster_IDHwt"]]),
    # list(cluster_IDHmut, snakemake@output[["cluster_IDHmut"]])
)



####################################################
####################################################


##### WRITE THE OUTPUTS

for (i in seq(1, length(clusters_results_list), 1)) {
    fwrite(
        x = as.data.table(clusters_results_list[[i]][1]),
        file = unlist(clusters_results_list[[i]][2]),
        sep = ","
    )
}









