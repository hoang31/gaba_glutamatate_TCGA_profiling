

#############################################
## GO and KEGG enrichment analysis of gene sets
#############################################

#############################################
## Load the libraries
#############################################

library(data.table)
library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(EnhancedVolcano)

#############################################
## Load the data
#############################################

## load the migration markers
migration_genes <- fread(
    snakemake@input[["migration_genes"]]
)

## load the deseq2 results of MIXED vs OTHERS
DE_res <- fread(
    snakemake@input[["DEgene_results_dir"]]
)

## set the cutoff for the pvalues
pval_cutoff <- 0.05
fc_cutoff <- 2

## extract only the genes that are differentially expressed between the two conditions
DE_res_filtered <- DE_res[padj < pval_cutoff,]
DE_res_filtered <- DE_res_filtered[(log2FoldChange < -fc_cutoff) | (log2FoldChange > fc_cutoff),]
DE_genes <- unlist(DE_res_filtered[, Gene_Name])


#############################################
## Enrichment analysis 
#############################################

## transform the gene symbol to other gene id
DE_genes <- bitr(
    DE_genes,
    fromType = "SYMBOL",
    toType = c("ENSEMBL", "ENTREZID"),
    OrgDb = 'org.Hs.eg.db'
)

DE_genes <- unlist(DE_genes[, 'ENTREZID'])

#####
## GO enrichment
#####

ego_CC <- enrichGO(
    gene = DE_genes,
    #universe = names(geneList),
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "bonferroni",
    pvalueCutoff = pval_cutoff,
    qvalueCutoff = 0.05,
    readable = TRUE
)

ego_MF <- enrichGO(
    gene = DE_genes,
    #universe = names(geneList),
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "bonferroni",
    pvalueCutoff = pval_cutoff,
    qvalueCutoff = 0.05,
    readable = TRUE
)

ego_BP <- enrichGO(
    gene = DE_genes,
    #universe = names(geneList),
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "bonferroni",
    pvalueCutoff = pval_cutoff,
    qvalueCutoff = 0.05,
    readable = TRUE
)

ego_CC <- as.data.table(summary(ego_CC))[order(Count, decreasing = T),]
ego_MF <- as.data.table(summary(ego_MF))[order(Count, decreasing = T),]
ego_BP <- as.data.table(summary(ego_BP))[order(Count, decreasing = T),]

#####
## KEGG enrichment
#####

kk <- enrichKEGG(
    gene = DE_genes,
    organism = 'hsa',
    pvalueCutoff = 0.05

)

kk <- as.data.table(summary(kk))[order(Count, decreasing = T),]


#############################################
## Comparison with the migration genes
#############################################


### extract the migrations genes that are differentially expressed
##DE_migration_genes <- migration_genes[gene_symbol %in% DE_genes]

#color_palette <- c(
#    brewer.pal(9, "Set1"),
#    brewer.pal(8, "Dark2")
#)

### merge the DESEQ2 information with the go terne annotation
#DE_gene_information <- merge(
#    DE_res,
#    migration_genes[, -c('go_id'), with = F],
#    by.x = 'Gene_Name',
#    by.y = 'gene_symbol',
#    sort = F,
#    all.x = T
#)[order(go_name)]

### add no migration information and color
##DE_gene_information[is.na(DE_gene_information[,go_name]), go_name := 'no_migration']

### put genes that are not associated with migration in grey
#DE_gene_information[is.na(DE_gene_information[,go_name]), color := 'grey']

### put genes that are not significant in grey
#DE_gene_information[
#    (abs(log2FoldChange) < fc_cutoff) | (padj > pval_cutoff),
#    color := 'grey'
#]

### associated a color with the go term
#color_table <- data.table(
#    go_name = names(table(DE_gene_information[is.na(DE_gene_information[,color]),go_name]))
#)
#color_table[, color := color_palette[1:(nrow(color_table))]]


#DE_gene_information[
#    is.na(DE_gene_information[,color]),
#    color := sapply(
#        go_name,
#        function(x) {
#            return(unlist(color_table[go_name == x,color]))
#        }
#    )
#]

#DE_gene_information[color == 'grey', go_name := 'ns']





#DE_gene_information <- DE_gene_information[ !(color == 'grey'),]
##DE_gene_information <- DE_gene_information[(go_name != 'neuron_migration'), color := 'grey']
##DE_gene_information <- DE_gene_information[(go_name != 'cell_migration'), color := 'grey']
##DE_gene_information <- DE_gene_information[(go_name != 'negative_regulation_of_cell_migration'), color := 'grey']
#DE_gene_information <- DE_gene_information[(go_name != 'positive_regulation_of_cell_migration'), color := 'grey']



#print(DE_gene_information)

### create the vector of color
#color_vector <- unlist(DE_gene_information[, color])
#names(color_vector) <- unlist(DE_gene_information[, go_name])

#table(names(color_vector))

#volcano_plot <- EnhancedVolcano(
#    toptable = DE_gene_information,
#    lab = unlist(DE_gene_information[, Gene_Name]),
#    selectLab = c(''),
#    x = "log2FoldChange",
#    y = "padj",
#    colCustom = color_vector,

#    title = "",
#    pCutoff = pval_cutoff,
#    FCcutoff = fc_cutoff,

#    colAlpha = 0.60,
#    pointSize = 8,
#    labSize = 5
#) + 
#    theme_minimal(base_size = 28) +
#    labs(color='Specie') +
#    labs(subtitle = 'Volcano Plot')

### save the graph
#ggsave(
#    'data/migration_volcanoplot.png',
#    plot = volcano_plot,
#    device = "png",
#    height = 10,
#    width = 20
#)

#Sys.sleep(100000)
#print(volcano_plot)
#break



##break


### merge this genes with the DESEQ2 informations such as fold change etc
#DE_migration_genes <- merge(
#    DE_migration_genes,
#    DE_res,
#    by.x = 'gene_symbol',
#    by.y = 'Gene_Name'
#)[order(go_name), ]

#print(DE_migration_genes)

##break




#############################################
## Write the files
#############################################


fwrite(
    ego_CC,
    'data/ego_CC.csv',
    sep = ','
)

fwrite(
    ego_MF,
    'data/ego_MF.csv',
    sep = ','
)

fwrite(
    ego_BP,
    'data/ego_BP.csv',
    sep = ','
)

fwrite(
    kk,
    'data/kk.csv',
    sep = ','
)

#break