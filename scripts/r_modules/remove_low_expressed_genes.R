
########################################################################
########################################################################

## Create heatmap from the expressiom data

########################################################################
########################################################################

## load packages
library(data.table)


########################################################################
########################################################################

##### load data

## load emsembl and symbols id associated with the GABA, GLUTA and CALCIUM pathways
gene_id <- fread(
    snakemake@params[["gene_id"]],
    sep = ",",
    header = T
)

## load kegg pathway genes
suppressMessages(genes <- scan(
    file = snakemake@input[["kegg_pathway_genes"]],
    what = character(),
    sep = "\n"
))

## load expression data
normalized_expression_data <- fread(
    file = snakemake@input[["normalized_expression_data"]],
    sep = ",",
    header = T
)

## load filter params
expression_filter <- as.numeric(snakemake@params[["expression_filter"]])

## load external scripts
source(snakemake@params[["utils_R"]])


########################################################################
########################################################################

##### PARSE THE EXPRESSION DATA

## rename the gene name column with "Gene_Name" and remove numbers after the "."
normalized_expression_data <- normalized_expression_data[
    , setnames(.SD, "rn", "Gene_Name")
    ][, Gene_Name := sapply(Gene_Name, function(x) strsplit(x, split = "[.]")[[1]][1])]

## replace the ensembl id to gene namess
normalized_expression_data <- normalized_expression_data[, Gene_Name := ensemblID_to_geneSymbole(unlist(normalized_expression_data[, "Gene_Name"]), gene_id)] 


########################################################################
########################################################################


##### FILTRATION STEPS

## extract the expression data associated with the GLUTA and GABA genes
normalized_expression_data_GABA_GLUTA <- normalized_expression_data[Gene_Name %in% genes, ]

## remove the genes whose the count max are lower than xxx
normalized_expression_data_expressionFiltered <- normalized_expression_data_GABA_GLUTA[
    apply(
        normalized_expression_data_GABA_GLUTA[, 2:ncol(normalized_expression_data_GABA_GLUTA)],
        1,
        function(x) max(as.numeric(x))
    ) >= expression_filter, 
]

## keep the genes whose 50% of the samples have a value superior than the expression_filter
#normalized_expression_data_expressionFiltered <- normalized_expression_data_expressionFiltered[
#    apply(
#        normalized_expression_data_expressionFiltered[, 2:ncol(normalized_expression_data_GABA_GLUTA)],
#        1,
#        function(x) sum(x > 5)
#    ) >= ncol(normalized_expression_data_expressionFiltered)*0.50, 
#]

########################################################################
########################################################################

##### RETURN OUTPUTS


## gene of interest that we retreive after the filtration of low expressed genes
genes_of_interest <- unlist(normalized_expression_data_expressionFiltered[, "Gene_Name"])

## extact the ensembl id of the interest genes
genes_of_interest <- gene_id[gene_name %in% genes_of_interest, ]

## return the genes of interest
fwrite(
    x = genes_of_interest,
    file = snakemake@output[["genes_highExpressed"]],
    sep = ",",
    col.names = T
)

## return the expression data
fwrite(
    x = normalized_expression_data_expressionFiltered,
    file = snakemake@output[["normalized_expression_data_expressionFiltered"]],
    sep = ",",
    col.names = T
)

## return the sample data
d_sample_data <- as.data.table(
    colnames(normalized_expression_data_expressionFiltered)[!(colnames(normalized_expression_data_expressionFiltered) == "Gene_Name")]
)

d_sample_data <- d_sample_data[, tissus := rep(
        "brain",
        nrow(d_sample_data)
),][,setnames(
        .SD,
        "V1",
        "samples"
)]

fwrite(
    x = d_sample_data,
    file = snakemake@output[["normalized_expression_data_sample_data"]],
    sep = ",",
    col.names = T
)



##########################################################################################################################################################################################################################################################