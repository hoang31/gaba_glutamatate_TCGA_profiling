
###################################################
###################################################


## Create heatmap from the expressiom data

###################################################
###################################################


## load packages
library(data.table)
suppressMessages(library(DESeq2))


print('okkkkk')
###################################################
###################################################


##### load data

## load expression data
# d_expression_CGGA = fread(file = snakemake@input[["CGGA_expression"]], sep = ",", header = T)
# d_expression_IGAP = fread(file = snakemake@input[["IGAP_expression"]], sep = ",", header = T)
d_expression_TCGA_cancer = fread(file = snakemake@input[["TCGA_expression"]], sep = ",", header = T)
d_expression_TCGA_normal = fread(file = snakemake@input[["TCGA_expression_normal"]], sep = ",", header = T)
# d_expression_GTEX = fread(file = snakemake@input[["GTEX_expression"]], sep = ",", header = T)

## load clinical data
# d_clinical_CGGA = fread(file = snakemake@input[["CGGA_clinical"]], sep = ",", header = T)
# d_clinical_IGAP = fread(file = snakemake@input[["IGAP_clinical"]], sep = ",", header = T)
d_clinical_TCGA_cancer = fread(file = snakemake@input[["TCGA_clinical"]], sep = ",", header = T)
d_clinical_TCGA_normal = fread(file = snakemake@input[["TCGA_clinical_normal"]], sep = ",", header = T)
d_clinical_TCGA_external = fread(snakemake@params[["TCGA_clinical_EXTERNAL"]], sep = ",", header = T)
# d_clinical_GTEX = fread(file = snakemake@input[["GTEX_clinical"]], sep = ",", header = T)

## laod exetrnal scripts
source(snakemake@params[["utils_R"]])


###################################################
###################################################


##### TCGA DATA


###################################################

##### Heatmap for all the glutamate and GABA genes

## rbind the TCGA data associated with cancer data and normal data
d_clinical_TCGA = rbind(d_clinical_TCGA_normal, d_clinical_TCGA_cancer)

## merge the expression data of cancer and normal data
d_expression_TCGA = merge(d_expression_TCGA_normal, d_expression_TCGA_cancer, by = "genes", all = T)
colnames(d_expression_TCGA)[1] = "Gene_Name" # rename first column

## remove the __alignment_not_unique, __ambiguous and __no_feature rows
row_to_remove = c("__alignment_not_unique", "__ambiguous", "__no_feature")
d_expression_TCGA = d_expression_TCGA[!(Gene_Name %in% row_to_remove), ]

## extract the id for the count files
d_clinical_TCGA = unique(d_clinical_TCGA[grep(x = d_clinical_TCGA[,File_Name], pattern = "counts"),])

## create a column called "type" associated with the tumor type, normal or cancer 
d_clinical_TCGA[,"type" := sapply(d_clinical_TCGA[,Sample_Type], function(x) ifelse(x == "Solid Tissue Normal", "normal", "cancer"))]

## remove duplicate column (Project ID)
d_clinical_TCGA = d_clinical_TCGA[,unique(names(d_clinical_TCGA)),with=FALSE]

## merge clinical data with other clinical
d_clinical_TCGA = merge(d_clinical_TCGA, 
                        d_clinical_TCGA_external, 
                        by.x = "Case_ID", 
                        by.y = "Case", 
                        all.x = T,
                        sort = F)

## create a column corresponding to the database name
d_clinical_TCGA[,"database" := "TCGA"]


###################################################
###################################################


##### Function for formate the data for the deseq2 object

## function for create the "coldata" data frame containing all the information of the samples
## example of EXPRESSION_MATRIX_INPUT :
## Gene_Name    sample1   sample2   sample3
## gene1        x         x         x
## gene2        x         x         x
## gene2        x         x         x

create_coldata = function(EXPRESSION_MATRIX_INPUT)
{
  coldata_output = data.frame(colnames(EXPRESSION_MATRIX_INPUT)[colnames(EXPRESSION_MATRIX_INPUT) != "Gene_Name"],
                       row.names = colnames(EXPRESSION_MATRIX_INPUT)[colnames(EXPRESSION_MATRIX_INPUT) != "Gene_Name"])
  colnames(coldata_output)[1] = "sample" # rename the column
  
  return(coldata_output)
}


## function for extracing the expression data related to the samples in the coldata
extract_expression_data = function(expression_matrix, col_data)
{
  expression_matrix_input = copy(expression_matrix)
  ## extract the exprssion related to the samples in the coldata
  expression_matrix_input = expression_matrix_input[,append("Gene_Name", rownames(col_data)), with = F]
  
  ## convert the element to numeric and then to data frame 
  setDF(expression_matrix_input)
  expression_output = data.frame(apply(expression_matrix_input[,2:ncol(expression_matrix_input)],c(1,2),
                                       function(x) round(x)), 
                                 row.names = expression_matrix_input[,"Gene_Name"],
                                 check.names = F)
  
  return(expression_output)
  
}


## function for apply DESEQ2 normalization 
deseq2_analysis = function(CTS_INPUT, COLDATA_INPUT)
{
  ## create the dds object with the cts and coldata data
  dds = DESeqDataSetFromMatrix(countData = CTS_INPUT,
                               colData = COLDATA_INPUT,
                               design = ~ 1)
  
  ## remove the genes not expressed
  keep = rowSums(counts(dds)) >= 10
  dds = dds[keep,]
  
  ## differential expression analysis
  dds = DESeq(dds, 
              parallel = F
              )
  
  ## estimate the size factor for the normalization
  dds <- estimateSizeFactors(dds)
  
  ## variance stabilizing transformation
  vsd = assay(varianceStabilizingTransformation(dds, blind = FALSE, fitType = "parametric"))
  
  ## extract the normalized counts
  normalized_cts = data.table(vsd, keep.rownames = T)
  
  return(normalized_cts)
}

###################################################
###################################################


##### Normalize the expression data by the function vst() available in the DESEQ data
# d_expression_TCGA = d_expression_TCGA[,1:10]

## create the "coldata", containing the informations of the samples
coldata = create_coldata(d_expression_TCGA)

## create the "cts" data frame containing the expression data which are associated with the id sample from the "coldata" data
cts = extract_expression_data(d_expression_TCGA, coldata)

## normalization of count data

TCGA_expression_normalized = deseq2_analysis(cts, coldata)
fwrite(TCGA_expression_normalized, snakemake@output[["TCGA_expression_normalized"]], sep = ",")


###################################################
###################################################