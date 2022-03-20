
###############################
## Filter rna proteins for reduce the file size
###############################

###############################
## load libraries
###############################

library(data.table)

###############################
## load the data
###############################

## load the interaction data
interaction_data <- fread(
    snakemake@input[["rnainter_rna_protein"]],
    sep = "\t"
)

## load the cutoff used for filtering
cutoff <- snakemake@params[["confidence_score_cutoff"]]

###############################
## Filter 
###############################

## set the references for the data table
setkey(
    x = interaction_data,
    "V4",
    "V5",
    "V10"
)

## filter by the organism
interaction_data <- interaction_data[V5 == "Homo sapiens",]

## filter by the molecule type
interaction_data <- interaction_data[V4 == "mRNA",]

## filter by the confidence score
interaction_data <- interaction_data[V10 > cutoff,]

###############################
## write the data
###############################

fwrite(
    interaction_data,
    snakemake@output[["rnainter_rna_protein_filtered"]],
    sep = "\t"
)