
###################################################
###################################################


##### COPY NUMBER VARIATION ANALYSIS


###################################################
###################################################


##### LOAD PACKAGES

library(data.table)
library(stringr)
library(magrittr)

###################################################
###################################################


##### LOAD INPUT DATA

dt_cnv_gbm <- fread(
    input = snakemake@input[["cnv_gbm"]],
    sep = "\t",
    header = T
)

## load the cnv data of LGG project
dt_cnv_lgg <- fread(
    input = snakemake@input[["cnv_lgg"]],
    sep = "\t",
    header = T
)

## load the cnv id
dt_cnv_id <- fread(
    input = snakemake@input[["cnv_id"]],
    sep = "\t",
    header = T
)

## load the gene_id
gene_id <- fread(
    input = snakemake@input[["gene_id"]],
    sep = ",",
    header = T
)

## load clinical data for extracting the sample id
dt_clinical_data_cancer <- fread(
    input = snakemake@input[["clinical_data"]],
    sep = ",",
    header = T
)

## remove the duplicates
dt_clinical_data_cancer <- dt_clinical_data_cancer[
        grep(
            x = dt_clinical_data_cancer[, File_Name],
            pattern = "counts"
        ),
    ]

## load "utils" functions
source(snakemake@params[["utils_R"]])


###################################################
###################################################


##### CNV ANALYSIS (copy number vatiation)


##### FUNCTIONS

## function for transform the ensembl id to gene symbol and remove gene that were not transformed
transform_and_extract_genes_cnv <- function(
    CNV_DATA_INPUT, # cnv data which posess a calumn called "Gene Symbol" corresponding to the ensembl id
    GENE_ID_INPUT # corresponding to a data table with a "gene_id" (ensembl ID) and "gene_name" (gene symbol) of interest genes
) {
    ## copy the input data
    cnv_data <- copy(CNV_DATA_INPUT)

    ## transforme the ensembl id associated with the gaba and gluta id (given by the "id" input) to gene symbols
    cnv_data[, "Gene Symbol" := ensemblID_to_geneSymbole(`Gene Symbol`, GENE_ID_INPUT)]

    ## extract only genes of interest
    cnv_data <- cnv_data[!grep(unlist(cnv_data[,"Gene Symbol"]), pattern = "ENS"),]
    
    ## return the  filtered cnv data
    return(cnv_data)
}


## function for counting the chromosomal aberrations of the data
## return a matrix of the counting
count_aberrations <- function(
    CNV_DATA_INPUT # data table containing in column all the sample names (+ "Gene Synbol" + "Gene ID" + "Cytoband") and the rows corresponding to genes
) {
    ## copy the cnv data input
    cnv_data <- copy(CNV_DATA_INPUT)

    ## extract the informations of each row
    information_data <- copy(cnv_data[, c("Gene Symbol", "Gene ID", "Cytoband"), with = F])

    ## extract the aberration informations for each samples
    cnv_data <- (cnv_data[, !c("Gene Symbol", "Gene ID", "Cytoband"), , with = F])

    ## initialize the counting data with the first counting
    counting_data <- as.data.table(
            table(unlist(cnv_data[1,]))
        )[, setnames(.SD, c("V1", "N"), c("aberration", as.character(information_data[1,"Gene Symbol"])) )]

    ## transform the aberration column to character
    counting_data <- counting_data[, aberration := as.character(aberration)]
    levels(counting_data[,"aberration"]) <- c("-1", "0", "1") # set the level

    ## for each line, do the counting of the aberrations
    for (i in seq(2, nrow(cnv_data), 1)) {

        ## do the counting
        counting <- as.data.table(
            table(unlist(cnv_data[i,]))
        )[, setnames(.SD, c("V1", "N"), c("aberration", as.character(information_data[i,"Gene Symbol"])))]

        ## transform the aberration column to character
        counting <- counting[, aberration := as.character(aberration)]
        levels(counting[,"aberration"]) <- c("-1", "0", "1")

        # print("------")
        # print(i)
        # print(as.character(information_data[i,"Gene Symbol"]))
        # print(counting)

        # merge the counting of the gene i with the counting data table
        counting_data <- merge(
            x = counting_data,
            y = counting,
            by = "aberration",
            sort = FALSE,
            all.x = TRUE
        )
    }

    ## transform to dataframe and set the rowname with tne aberration column
    counting_data <- data.frame(
        counting_data, 
        row.names = counting_data[,aberration]
    )

    ## remove the aberration column
    counting_data[, "aberration"] <- NULL

    ## transpose the matrix
    counting_data <- t(counting_data)
    
    ## transform to data table
    counting_data <- as.data.table(counting_data, keep.rownames = T)

    ## transform NA values to 0
    counting_data[is.na(counting_data)] <- 0
    
    ## return the counting data
    return(counting_data)

}







##### CNV ANALYSIS (copy number vatiation)

## remove the gene ID and Cytoband columns
dt_cnv_lgg[, c("Gene ID", "Cytoband") := NULL]

## merge the GBM and LGG cnv data
dt_cnv_all <- merge(
    x = dt_cnv_gbm,
    y = dt_cnv_lgg,
    by = "Gene Symbol",
    all.x = TRUE,
    sort = FALSE
)

## extract only the file name without extension
dt_cnv_id[, "File Name" := sapply(
    `File Name`,
    function(x) strsplit(
        x,
        split = "[.]"
    )[[1]][2]
)]

## remove case ID and Sample ID duplicates
dt_cnv_id[, "Case ID" := sapply(
    `Case ID`,
    function(x) strsplit(
        x,
        split = "[,]"
    )[[1]][1]
)][, "Sample ID" := sapply(
    `Sample ID`,
    function(x) strsplit(
        x,
        split = "[,]"
    )[[1]][1]
)]

## take only the File Name and Case ID
dt_cnv_id <- dt_cnv_id[, c("File Name", "Case ID")]

## extract the cnv id associated with the samples that possess rna seq data
dt_cnv_id2 <- dt_cnv_id[`Case ID` %in% unlist(dt_clinical_data_cancer[, c("Case_ID"), with = F]),]

## extract samples possessing several duplicates
samples_with_replicates <- data.table(table(dt_cnv_id2[, "Case ID"]))[N > 1,V1]

## extract the sample replicates positions
for (i in samples_with_replicates) {
    sample_id <- grep(
        x = unlist(dt_cnv_id2[, "Case ID"]),
        pattern = i
    )
    sample_id <- sample_id[2:length(sample_id)]
    dt_cnv_id2[sample_id, "File Name" := "toremove"]
}

## remove the replicates
dt_cnv_id2 <- dt_cnv_id2[`File Name` != "toremove",]

## extract the cnv data if the sample which possess rna-seq data
dt_cnv_all2 <- dt_cnv_all[,append(c("Gene Symbol", "Gene ID", "Cytoband"), unlist(dt_cnv_id2[, "File Name"])), with = F]

## rename the column (corresponding to the file name) by the id sample
dt_cnv_all2 <- dt_cnv_all2[, setnames(
    .SD,
    colnames(dt_cnv_all2),
    sapply(
        colnames(dt_cnv_all2),
        function(x) {
            id_sample <- unlist((dt_cnv_id2[`File Name` == x, "Case ID"]))
            names(id_sample) <- NULL
            if (length(id_sample) == 0) {return(x)}
            else {return(id_sample)}
        }
    )
)
]



#############################################
## Extract alteration of interest
#############################################

## id of gene of interest
gene_of_interest <- c(
    'EGFR' = 'ENSG00000146648',
    'CDKN2A' = 'ENSG00000147889',
    'CDKN2B' = 'ENSG00000147883',
    'MYB' = 'ENSG00000118513',
    'MYBL1' = 'ENSG00000185697',
    'FGFR' = 'ENSG00000077782'
)

## collapse all the id together
gene_of_interest_collapsed <- paste(
    gene_of_interest,
    collapse = '|',
    sep = ''
)

## extract the data related to the genes of interest
alterations_of_interest <- dt_cnv_all2[grep(`Gene Symbol`, pattern = gene_of_interest_collapsed), ]

## change the ensembl id to gene symboles
alterations_of_interest <- alterations_of_interest[
    ,
    `Gene Symbol` := sapply(
        `Gene Symbol`,
        function(x) {

            ## remove the number after the point
            splitted_id <- str_split(
                string = x,
                pattern = '[.]'
            )[[1]][1]
            
            ## retreive the names associated with the ensembl id
            gene_symb <- names(
                    gene_of_interest[
                        grep(
                        gene_of_interest,
                        pattern = splitted_id
                    )
                ]
            )

            return(gene_symb)
        }
    )
]


#############################################
#############################################


## extract the genes of interest and transform the emsembl id to gene symbols
dt_cnv_all2 <- transform_and_extract_genes_cnv(dt_cnv_all2, gene_id)

dt_cnv_all2 <- rbind(
    dt_cnv_all2,
    alterations_of_interest
)
## do the counting of all the chromosomal aberrations
CNV_counting <- count_aberrations(dt_cnv_all2)

## transform the -1 and 1 to "deletion" and "amplification"
TCGA_cnv_data <- copy(dt_cnv_all2)
TCGA_cnv_data[TCGA_cnv_data == 1] <- "amplification"
TCGA_cnv_data[TCGA_cnv_data == -1] <- "deletion"
TCGA_cnv_data[TCGA_cnv_data == 0] <- "0"

## extract the cnv information
alterations_of_interest <- copy(TCGA_cnv_data)
alterations_of_interest <- alterations_of_interest[
    ,
    !c('Gene ID', 'Cytoband'),
    with = F
]

## do the transposition of the alteration data table
alterations_of_interest <- transpose(
    rbindlist(
        list(
            as.list(colnames(alterations_of_interest)),
            alterations_of_interest
        )
    )
)

## rename the colnames by the alteration
alterations_of_interest <- alterations_of_interest[
    ,
    setnames(
        .SD,
        colnames(alterations_of_interest),
        unlist(alterations_of_interest[V1 == 'Gene Symbol',])
    )
][
    ,
    setnames(
        .SD,
        'Gene Symbol',
        'sample_id'
    )
][
    sample_id != 'Gene Symbol',
]

## return the data
fwrite(
    alterations_of_interest,
    snakemake@output[["alteration_data"]],
    sep = ",")











###################################################
###################################################


# ##### SNV ANALYSIS (single nucleotide variation)


# ##### FUNCTIONS

# ## function for extract the snv data of the genes of interest
# transform_and_extract_genes_snv <- function(
#     SNV_DATA_INPUT, # snv data which posess a calumn called "Hugo_Symbol" corresponding to the Hugo Symbol id
#     GENE_ID_INPUT # corresponding to a data table with a "gene_id" (ensembl ID) and "gene_name" (gene symbol) of interest genes
# ) {
#     ## copy the input data
#     snv_data <- copy(SNV_DATA_INPUT)

#     ## set keys
#     setkey(snv_data, Hugo_Symbol)

#     ## extract only genes of interest
#     snv_data <- snv_data[Hugo_Symbol %in% unlist(GENE_ID_INPUT[,"gene_name"]),]

#     ## return the  filtered cnv data
#     return(snv_data)
# }


# ## function for reading all the snv data contained in the directory
# read_all_snv_files <- function(
#     DIR_FILE_PATH_INPUT, # path of the directory containing the snc data
#      GENE_ID_INPUT # corresponding to a data table with a "gene_id" (ensembl ID) and "gene_name" (gene symbol) of interest genes
# ) {

#     ## read all the file from the dir file path input
#     snv_files <- list.files(snakemake@input[["snv_data"]])
    
#     ## create the path to acess to the snv files
#     snv_files_paths <- sapply(
#         snv_files,
#         function(x) paste(
#             DIR_FILE_PATH_INPUT,
#             x,
#             sep = ""
#         )
#     )

#     ## initialize the cnv data reading first file
#     snv_data <- fread(
#             input = snv_files_paths[1],
#             sep = "\t",
#             header = T
#         )
    
#     ## string of all the snv analysis methods
#     pattern_string <- "muse|mutect|somaticsniper|varscan"

#     ## put the name of the snv analysis methods
#     snv_data[, snv_methods := str_extract(
#         string = snv_files_paths[1],
#         pattern = pattern_string
#         )
#     ]

#     ## extract snv data only for the gene of interest (stocked in the "GENE_ID_INPUT" variable)
#     snv_data <- transform_and_extract_genes_snv(
#         SNV_DATA_INPUT = snv_data,
#         GENE_ID_INPUT = GENE_ID_INPUT
#     )

#     ## read the other files
#     for (i in seq(2,length(snv_files_paths), 1)) {
        
#         ## read the i-th file
#         snv <- fread(
#             input = snv_files_paths[i],
#             sep = "\t",
#             header = T
#         )

#         ## extract snv data of the genes of interest
#         snv <- transform_and_extract_genes_snv(
#             SNV_DATA_INPUT = snv,
#             GENE_ID_INPUT = GENE_ID_INPUT
#         )

#         ## put the name of the snv analysis methods
#         snv[, snv_methods := str_extract(
#             string = snv_files_paths[i],
#             pattern = pattern_string
#             )
#         ]

#         ## rbind the snv data to the snv_data datatable
#         snv_data <- rbind(
#             snv_data,
#             snv
#         )
#     }

#     ## coloumn to keep from the results
#     col_to_keep <- c(
#     "snv_methods",
#     "Hugo_Symbol",
#     "Chromosome",
#     "Start_Position",
#     "End_Position",
#     "Variant_Classification",
#     "Consequence",
#     "Tumor_Sample_Barcode"
#     )

#     ## extract the columns of interest
#     snv_data <- snv_data[,col_to_keep, with = F]

#     snv_data[, Tumor_Sample_Barcode := sapply(
#         Tumor_Sample_Barcode,
#         function(x) str_sub(
#             x,
#             start = 1,
#             end = 12
#         )
#     )]

#     return(snv_data)
# }


# ##### SNV ANALYS

# ## get all the mutations occuring in the genes
# d_mutations <- read_all_snv_files(
#     DIR_FILE_PATH_INPUT = snakemake@input[["snv_data"]],
#     GENE_ID_INPUT = gene_id
# )

# ## extract the mutation associated with a least two snv methods
# d_mutations_filtered <- d_mutations[, .N, by = c("Tumor_Sample_Barcode", "Hugo_Symbol", "Start_Position", "End_Position")][N > 1, ]

# ## copy the structure of the expression data (we want as columns the samples and in rows the gene names)
# TCGA_cnv_snv_data <- copy(TCGA_cnv_data[, !c("Cytoband", "Gene ID")])





# ## extact all the mutation from the mutation data and initialize the data table containing all these informations
# for (i in seq(1, nrow(d_mutations_filtered), 1)) {
    
#     ## extract the gene name affected by the alteration (corresponding to the rows of the "TCGA_cnv_snv_data")
#     gene_name <- toString(d_mutations_filtered[i, "Hugo_Symbol"])
 
#     ## extract the sample id affected by the alteration (corresponding to the rows of the "TCGA_cnv_snv_data")
#     sample <- toString(d_mutations_filtered[i, "Tumor_Sample_Barcode"])

#     ## extract the position of the gene name    
#     gene_position <- grep(
#         x = unlist(TCGA_cnv_snv_data[,"Gene Symbol"]),
#         pattern = gene_name
#     )


#     # cat("---", gene_name, ",", sample, "---", "\n")



#     ## if the sample exist
#     if (sample %in% colnames(TCGA_cnv_snv_data) == TRUE & gene_name %in% unlist(TCGA_cnv_snv_data[, "Gene Symbol"])) {

#         ## if the sample associated with the gene possesses a deletion or a amplification, add "mutation" to the "deletion/amplification" string
#         if ((TCGA_cnv_snv_data[gene_position, sample, with = F] == "deletion") | (TCGA_cnv_snv_data[gene_position, sample, with = F] == "amplification")) {

            
#             # cat("---", gene_position, ",", sample, "---", "\n")
#             # set(
#             # x = TCGA_cnv_snv_data,
#             # i = gene_position,
#             # j = sample,
#             # value = paste("mutation", toString(TCGA_cnv_snv_data[gene_position, sample, with = F]), sep = "_")
#             # )

#             set(
#             x = TCGA_cnv_snv_data,
#             i = gene_position,
#             j = sample,
#             value = toString(TCGA_cnv_snv_data[gene_position, sample, with = F])
#             )
#         }
#         ## if not, just put "mutation"
#         else {
#             set(
#             x = TCGA_cnv_snv_data,
#             i = gene_position,
#             j = sample,
#             value = "mutation"
#             )
#         }
#     }
# }


# ###################################################
# ###################################################

# ##### ONCOPRINT

# ##### load the DE genes 

# DEgenes_directory <- snakemake@input[["DEgene_results_dir"]]

# ## retreive the file names
# file_names <- list.files(DEgenes_directory)

# ## retreive the file paths
# file_path <- sapply(
#     X = list.files(DEgenes_directory),
#     function(x) paste(
#         DEgenes_directory,
#         x,
#         sep = "/"
#     )
# )

# ## read all the files using the path files
# merged_file <- lapply(
#     X = file_path,
#     function(x) fread(
#         input = x,
#         sep = ","
#     )
# )

# ## transform the list to data table
# merged_file <- as.data.table(merged_file[[1]])

# ## initialize the cutoff for the fold change and the pvalue
# pval_cutoff <- 0.001
# FCcutoff <- 1

# ## extract genes and their informations which are below the pvalue cutoff
# # dt_DE_genes <- copy(merged_file[((log2FoldChange >= -FCcutoff) & (log2FoldChange <= FCcutoff)) | padj < pval_cutoff, "Gene_Name"])
# dt_DE_genes <- copy(merged_file[padj < pval_cutoff, "Gene_Name"])

# ## extract the DE genes
# DE_genes <- unlist(dt_DE_genes[, "Gene_Name"])

# ## extract the cnv and snv data of the DE genes
# TCGA_cnv_snv_data <- TCGA_cnv_snv_data[`Gene Symbol` %in% DE_genes,]




# d_cluster <- fread(
#     "~/Documents/gaba_and_glutamate_pathways_in_glioma/data/others/TCGA_IDHall_clusters.csv",
#     sep = ","
# )

# cluster_MIXTE_id <- d_cluster[cluster %in% c("MIXTE", "IDH_MUT_codel", "IDH_MUT_noncodel"),sample_id]



# ##### create the oncoprint
# library("ComplexHeatmap")

# d_alteraton_data <-  copy(TCGA_cnv_snv_data)


# d_alteraton_data <- d_alteraton_data[, colnames(d_alteraton_data) %in% append("Gene Symbol",cluster_MIXTE_id ), with = F]
# print(dim(d_alteraton_data))









# setDF(
#     d_alteraton_data,
#     rownames = unlist(d_alteraton_data[,"Gene Symbol"])
# )

# ## remove the Gene Symbol column
# d_alteraton_data[["Gene Symbol"]] <- NULL

# ## replace "0" value bye a blank
# d_alteraton_data[d_alteraton_data == "0"] <- ""

# ## colors for the oncoprint
# col <- c(
#     "deletion" = "blue3",
#     "amplification" = "red3",
#     "mutation" = "#008000",
#     "mutation_deletion" = "#7777fc",
#     "mutation_amplification" = "#ff6767"
# )

# col <- c(
#     "deletion" = "blue3",
#     "amplification" = "red3",
#     "mutation" = "#008000",
#     "mutation_deletion" = "blue3",
#     "mutation_amplification" = "red3"
# )

# alter_fun <- list(
#     background = function(x, y, w, h) {
#         grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
#             gp = gpar(fill = "#CCCCCC", col = NA))
#     },
#     # big blue
#     deletion = function(x, y, w, h) {
#         grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
#             gp = gpar(fill = col["deletion"], col = NA))
#     },
#     # big red
#     amplification = function(x, y, w, h) {
#         grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
#             gp = gpar(fill = col["amplification"], col = NA))
#     },
#     # small green
#     mutation = function(x, y, w, h) {
#         grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
#             gp = gpar(fill = col["mutation"], col = NA))
#     },

#     # 
#     mutation_deletion = function(x, y, w, h) {
#         grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
#             gp = gpar(fill = col["deletion"], col = NA))
#     },

#     # small green
#     mutation_amplification = function(x, y, w, h) {
#         grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
#             gp = gpar(fill = col["deletion"], col = NA))
#     }
# )


# ## create the oncoprint
# oncoprint1 <- oncoPrint(
#     d_alteraton_data,
#     alter_fun = alter_fun,
#     col = col,
#     column_title = "OncoPrint for TCGA GBM and LGG projects of GABA and Glutamate genes",
#     heatmap_legend_param = list(
#         title = "Alterations",
#         at = c("deletion", "amplification", "mutation"),
#         labels = c("Deletion", "Amplification", "Mutation")),
#     remove_empty_columns = TRUE,
#     remove_empty_rows = TRUE
# )

# jpeg(
#     filename =  "~/Documents/gaba_and_glutamate_pathways_in_glioma/data/figures/oncoprint.jpeg",
#     width = 800,
#     height = 1000
# )
# print(oncoprint1)
# dev.off()





# ###################################################
# ###################################################

# ##### OUTPUTS

# ## return the alteraton data table
# fwrite(
#     TCGA_cnv_snv_data,
#     snakemake@output[["alteration_data"]],
#     sep = ",")