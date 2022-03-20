
##########################################
## Calculate the correlation between IDH gene
##########################################

##########################################
## load the libraries
##########################################

library(data.table)
library(tidyverse)
source(snakemake@params[["utils_R"]])


##########################################
## load the data
##########################################


## load the expression data 
expression_data_dt <- fread(
    snakemake@input[["TCGA_expression_normalized"]]
)[, setnames(.SD, "rn", "gene_id")][, gene_id := lapply(gene_id, function(x) str_split(x, pattern = '[.]')[[1]][1])]

## load the cluster group data
cluste_group_dt <- fread(
    snakemake@input[["cluster_group"]],
    sep = ","
)



##########################################
## ANALYSIS
##########################################

## id for the IDH gene
id_IDH <- "ENSG00000138413"

## id for the genes of interest
id_of_interest <- c(
    #"ENSG00000138413", # IDH1
    #"ENSG00000182054", # IDH2
    #"ENSG00000169083",
    "ENSG00000146648"
)


## extract the genes of interest from the expression data
expression_data_filtered <- expression_data_dt[gene_id %in% id_of_interest, ]



## do the transposition of the expression_data_filtered
expression_data_filtered <- transpose_datatable(
    expression_data_filtered,
    column_name = "gene_id",
    new_name = 'sample_id'
)


#####
## Calcultate the correlation between the IDH1 and Androgen receptor gene
#####




## generate the point cloud for the visualization of the correlation
#correlation_plot <- ggplot(
#    data = expression_data_filtered,
#    aes(
#        x = as.numeric(ENSG00000138413),
#        y = as.numeric(ENSG00000169083)
#    )
#) + geom_point()

#print(correlation_plot)


## formate the data
expression_data_filtered <- melt(
    expression_data_filtered,
    id.vars = 'sample_id',
    measure.vars = id_of_interest
)

## merge with the cluster group data
expression_data_filtered <- merge(
    expression_data_filtered,
    cluste_group_dt,
    by = 'sample_id',
    all.y = T
)

#expression_data_filtered$cluster <- factor(
#    expression_data_filtered$cluster,
#    levels = c(
#        'MIXED',
#        'IDHmut_CODEL',
#        'IDHmut_nonCODEL',
#        'IDHwt'
#    ),

#)


## generate expression boxplot
expression_boxplot <- ggplot(
    expression_data_filtered,
    aes(
        x = cluster,
        y = as.numeric(value),
        fill = variable
    )
) +
    geom_boxplot(size = 3)
    #geom_bar(stat="identity", position=position_dodge())

print(expression_boxplot)


break



#####
## Normal samples
#####

### load the expression of controle samples
#controle_expression_dt <- fread(
#    "/home/hoang/Downloads/CGGA_RNAseq_Control_20.txt/CGGA_RNAseq_Control_20.txt",
#    sep = "\t"
#)

#genes_of_interest_vector <- c(
#    "IDH1",
#    "AR"
#)

### extract the expression data related to the genes of interest
#controle_expression_dt <- controle_expression_dt[gene_name %in% genes_of_interest_vector,]

### do the transposition
#controle_expression_dt <- transpose_datatable(
#    controle_expression_dt,
#    column_name = "gene_name",
#    new_name = 'sample_id'
#)

### formate the data
#controle_expression_dt <- melt(
#    controle_expression_dt,
#    id.vars = 'sample_id',
#    measure.vars = genes_of_interest_vector
#)

### add the cluster column 
#controle_expression_dt <- controle_expression_dt[, cluster := 'controle']

### rbind
##controle_expression_dt <- rbind(
##    controle_expression_dt,
##    expression_data_filtered
##)


#expression_boxplot <- ggplot(
#    controle_expression_dt,
#    aes(
#        x = variable,
#        y = as.numeric(value),
#        fill = cluster
#    )
#) +
#    geom_boxplot(size = 3)
#    #geom_bar(stat="identity", position=position_dodge())

#print(expression_boxplot)















break












test1 <- unlist(expression_data_dt[gene_id %in% id_IDH, !c('gene_id'), with = F])
test2 <- unlist(expression_data_dt[gene_id %in% id_of_interest, !c('gene_id'), with = F])


cor(
    test1,
    test2
)

test1 <- data.table(
    gene1 = test1,
    gene2 = test2
)



p <- ggplot(
    data = test1,
    aes(
        x = gene1,
        y = gene2
    )
) + 
    geom_point() +
    geom_smooth(method=lm)


print(p)




#print(expression_data_dt[1:5,1:3])
break