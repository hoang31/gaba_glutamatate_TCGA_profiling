
###################################################
## Extract the genes associated with the migrations
###################################################

###################################################
## Load the libraries 
###################################################

library(data.table)

###################################################
## Load the data
###################################################

## load the gene ontology id of interest associated with the migration
migration_go_id <- fread(
    snakemake@input[["migration_go_id"]],
    sep = ","
)

## load the gene ontology annotation associated with the homo sapiens specie
go_annotation <- fread(
    snakemake@input[["go_annotation"]],
    sep = "\t"
)

###################################################
## Analysis
###################################################

## rename the colnames of the go_annotation
colnames(go_annotation) <- c(
    "DB",
    "DB_Object_ID",
    "gene_symbol",
    "Qualifier",
    "go_id",
    "DB:Reference",
    "Evidence_Code",
    "With_(or)_From",
    "Aspect",
    "DB_Object_Name",
    "DB_Object_Synonym_(|Synonym)",
    "DB_Object_Type",
    "Taxon(|taxon)",
    "Date",
    "Assigned_By",
    "Annotation_Extension",
    "Gene_Product_Form_ID"
)

## extract the migration go id
go_id <- unlist(migration_go_id[, go_id])

## create the pattern to retreive all the go_id of interest 
go_id_regex <- paste(
    go_id,
    collapse = '|'
)

## retreive the genes that are associated with the go_id of interest
migration_dt <- go_annotation[grep(go_id, pattern = go_id_regex),]

## remove duplicated genes
migration_dt <- unique(migration_dt[, gene_symbol, by = 'go_id'])

## merge the migration data with the go_id
migration_dt <- merge(
    migration_dt,
    migration_go_id,
    by = 'go_id',
    all.x = T,
    sort = F
)

## for each gene, collapse the go id
migration_dt <- migration_dt[
    ,
    .(
        go_name = paste(go_name, collapse = ';'),
        go_id = paste(go_id, collapse = ';')

    ), 
    by = 'gene_symbol'
]


###################################################
## Write the outputs
###################################################

## write the migration marker data table
fwrite(
    migration_dt,
    snakemake@output[["migration_genes"]],
    sep = ","
)