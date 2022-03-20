

#######################################
## Generate a figure for the visualization of the alterations
#######################################

#######################################
## Load the libraries
#######################################

library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)

#######################################
## Load the data
#######################################

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

## load the epidemiological data
dt_epidemiological_data <- fread(
    snakemake@input[["epidemiological_data"]],
    sep = ',',
)

## load external scripts
source(snakemake@params[["heatmap_function"]])
source(snakemake@params[["clustering_metrics"]])
source(snakemake@params[["utils_R"]])

#######################################
## Create the heatmap
#######################################

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

## create a variable called "classification" describing the glioma classification : IDH-mut codel, IDH-mut noncodel and IDHwt
d_clinical_TCGA[`IDH status` == "WT", "classification" := "IDHwt"]
d_clinical_TCGA[(`IDH status` == "Mutant") & (`1p/19q codeletion` == "codel"), "classification" := "IDHmut_CODEL"]
d_clinical_TCGA[(`IDH status` == "Mutant") & (`1p/19q codeletion` == "non-codel"), "classification" := "IDHmut_nonCODEL"]


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

####################################################
####################################################
##### SEPARATION OF IDH-WT, IDH-MUT with 1p/19q codel and IDH-mut without the 1p/19q codel

## extract the sample assoiated with the glioma samples
d_clinical_cancer <- copy(d_clinical_TCGA[(type == "cancer") & classification %in% c("IDHwt", "IDHmut_CODEL", "IDHmut_nonCODEL"), ])

## extract the expresssion data associated with the sample of interest
d_expression_cancer <- copy(d_expression[
    ,
    append(
        "Gene_Name",
        unlist(d_clinical_cancer[, "Case_ID"])
    ),
    with = F
])


###################################################
##### GENERATE THE HEATMAPS 
###################################################

## extract the names of the duplicated columns
duplicated_column_names <- intersect(
    colnames(d_clinical_cancer),
    colnames(dt_epidemiological_data)
)
duplicated_column_names <- duplicated_column_names[duplicated_column_names != 'Case_ID']

## remove duplicated columns in the epidemiolocal data
set(dt_epidemiological_data, , duplicated_column_names, NULL)

## merge the clinical data with the epidemiological_data
d_clinical_cancer <- merge(
    d_clinical_cancer,
    dt_epidemiological_data,
    by = "Case_ID",
    all.x = TRUE,
    sort = F
)



## initialize the variables
var <- c(
    "site_of_resection_or_biopsy",
    "Transcriptome Subtype"


    #"Histology",
    #"Grade",

    #"IDH status",
    #"1p/19q codeletion",
    
    #"ATRX status",

    #"EGFR",
    #"TERT promoter status",
    #"Chr 7 gain/Chr 10 loss",

    #"CDKN2A",
    #"CDKN2B",

    #"cancer type"
)

## extract the clinical data related to the variable of intest
d_clinical_cancer <- d_clinical_cancer[, var, with = F]

## formate some legend names
d_clinical_cancer[d_clinical_cancer == "wt"] <- "WT"
d_clinical_cancer[d_clinical_cancer == "amplification"] <- "Amplification"
d_clinical_cancer[d_clinical_cancer == "deletion"] <- "Deletion"
d_clinical_cancer[d_clinical_cancer == "codel"] <- "Codel"
d_clinical_cancer[d_clinical_cancer == "non-codel"] <- "Non-codel"

d_clinical_cancer[d_clinical_cancer == "astrocytoma"] <- "Astrocytoma"
d_clinical_cancer[d_clinical_cancer == "glioblastoma"] <- "Glioblastoma"
d_clinical_cancer[d_clinical_cancer == "oligoastrocytoma"] <- "Oligoastrocytoma"
d_clinical_cancer[d_clinical_cancer == "oligodendroglioma"] <- "Oligodendroglioma"



## generate a heatmap of the cancer data generating the groups of interest
heatmap_IDHall <- create_heatmap_annotation(
    expression_data = d_expression_cancer,
    clinical = d_clinical_cancer,
    variables = var,
    #row_variables = gene_annotation,
    genes = unlist(d_expression_cancer[, Gene_Name]),
    pseudo_count = 0.1,
    path_file = snakemake@output[["alteration_visualization"]],
    log_normalization = F,
    col_split = 4,
    row_split = 3,
    clustering_distance_rows = "minkowski",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "spearman",
    clustering_method_columns = "ward.D2",
    horizontal_legend = F
)

## generate a heatmap of the cancer data generating the groups of interest
heatmap_legend <- create_heatmap_annotation(
    expression_data = d_expression_cancer,
    clinical = d_clinical_cancer,
    variables = var,
    #row_variables = gene_annotation,
    genes = unlist(d_expression_cancer[, Gene_Name]),
    pseudo_count = 0.1,
    path_file = snakemake@output[["alteration_visualization_legend"]],
    log_normalization = F,
    col_split = 4,
    row_split = 3,
    clustering_distance_rows = "minkowski",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "spearman",
    clustering_method_columns = "ward.D2",
    horizontal_legend = T
)


#print("finiiii")
#Sys.sleep(1000000000)