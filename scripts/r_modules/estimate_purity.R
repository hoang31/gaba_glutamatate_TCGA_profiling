

##########################################
##########################################


##### install R packages that are not available with conda


##########################################
##########################################


##### LOAD LIBRARIES
library(data.table)
library(estimate)
library(tidyverse)

##########################################
##########################################


##### LOAD DATA

## load util function
source(snakemake@params[["utils_R"]])

## load the ensembl id
ensembl_id <- fread(
    input = snakemake@params[["ensembl_id"]],
    sep = ","
)

# load FPKM expression data
fpkm_expression_data <- fread(
    input = snakemake@input[["fpkm_expression_data"]],
    sep = ","
)

# load FPKM data of the five normal tissues
normal_fpkm_expression_data <- fread(
    input = snakemake@input[["normal_fpkm_expression_data"]],
    sep = ","
)


##########################################
## transformation of sata
##########################################


## function for calculating the tumor purity of samples, taking as input the FPKM count matrix
calculate_purity <- function(
    fpkm_count_matrix, # the fpkm matrix that contains in row the genes and in column the samples. The gene column has to be called "genes"
    file_path # path where the file will be saved 
) {
    
    ## copy the fpkm data into a new data table for avoiding the dependancies
    fpkm_count_matrix_input <- copy(fpkm_count_matrix)
    
    ## tranform the ensembl id to gene symbol and extract the gene of interest from the gene_id data
    fpkm_count_matrix_input[, genes := ensemblID_to_geneSymbole(genes, ensembl_id)]
    fpkm_count_matrix_input <- fpkm_count_matrix_input[!grep(genes, pattern = "ENS"),][, setnames(.SD, "genes", "Gene_Name")]

    ## find the replicated genes
    dupplicated_genes <- unlist(as.data.table(table(unlist(fpkm_count_matrix_input[,"Gene_Name", with = F])))[order(N),][N > 1,V1])

    ## remove the duplicated genes from the fpkm expession data
    fpkm_count_matrix_input <- fpkm_count_matrix_input[!(Gene_Name %in% dupplicated_genes),]

    ## transform the expression data table to data frame and rename the row name with the gene column and then suppress the "Gene_Name" column
    setDF(fpkm_count_matrix_input, rownames = unlist(fpkm_count_matrix_input[,"Gene_Name", with = F]))
    fpkm_count_matrix_input[["Gene_Name"]] <- NULL

    ## write the expression file 
    fwrite(
        x = fpkm_count_matrix_input,
        file = "purity_expression.tsv",
        sep = "\t",
        row.names = T
    )

    ## extract the genes of interest
    filterCommonGenes(
        input.f= "purity_expression.tsv",
        output.f="purity.csv",
        id="GeneSymbol"
    )

    # estimate the purity of each sample
    estimateScore(
        "purity.csv",
        "purity.gct",
        platform = "affymetrix" 

    )

    ## read the purity estimation scores
    estimation_score_data <- fread(
        input = "purity.gct",
        sep = "\t",
        skip = 2
    )

    ## rename the sample replacing the '.' by "-s"
    colnames(estimation_score_data) <- str_replace_all(
        string = colnames(estimation_score_data),
        pattern = "[.]",
        replacement = "-"
    )

    ## add the column name as the first row and do the transposition of the matrix
    estimation_score_data <- as.data.table(
        transpose(
            rbindlist(
                list(
                    as.data.table(as.list(colnames(estimation_score_data))),
                    estimation_score_data
                )
            )
        )
    )

    ## rename the colnames and remove two row associated with the sample id
    estimation_score_data <- estimation_score_data[
        ,
        setnames(
            .SD,
            colnames(estimation_score_data),
            unlist(estimation_score_data[V1 == "NAME", ])
        )
    ][
        !(NAME == "NAME"),
    ][
        !(NAME == "Description"),
    ]

    ## delete the intermediate files
    unlink("purity_expression.tsv")
    unlink("purity.csv")
    unlink("purity.gct")

    ## write the datatable containing the tumor purity estimation score data
    fwrite(
        x = estimation_score_data,
        file = file_path,
        sep = ","
    )
}


calculate_purity(
    fpkm_expression_data,
    snakemake@output[["purity_score_tumor"]]
)


calculate_purity(
    normal_fpkm_expression_data,
    snakemake@output[["purity_score_normal"]]
)





###############################################
###############################################
## test of the purity calculation on the healthy tissus of the CGGA

#CGGA_healthy_samples <- fread(
#    '/home/hoang/Documents/PROJECTS/gaba_and_glutamate_pathways_in_glioma_V2/data/CGGA_RNAseq_Control_20.txt/CGGA_RNAseq_Control_20.txt',
#    sep = "\t"
#)

### rename the gene name column
#CGGA_healthy_samples <- CGGA_healthy_samples[
#    ,
#    setnames(
#        .SD,
#        c('gene_name'),
#        c('Gene_Name')
#    )
#]


#calculate_purity(
#    CGGA_healthy_samples,
#    '/home/hoang/Documents/PROJECTS/gaba_and_glutamate_pathways_in_glioma_V2/data/CGGA_RNAseq_Control_20.txt/CGGA_tumor_purity_normal.csv'
#)

#print(dim(CGGA_healthy_samples))
#print(colnames(CGGA_healthy_samples))
#print(CGGA_healthy_samples[1:3, 1:3])


###############################################
###############################################





#print('finiiiiiiiii')
#Sys.sleep(1000)