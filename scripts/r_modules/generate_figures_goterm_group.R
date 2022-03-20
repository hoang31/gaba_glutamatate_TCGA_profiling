
########################################
## Generate subgroup of GO terms together
########################################

########################################
## Load the libraries
########################################


options(width=system("tput cols", intern=TRUE))
library(data.table)
library(tidyverse)
library(ggalluvial)
library(egg)


########################################
## Load the data
########################################


## load the go group data 
go_groups_dt <- fread(
    snakemake@input[['go_groups']],
    sep = ','
)

## load some utils fonction
source(snakemake@params[['utils_R']])

## create the directory that will contain the data
dir.create(snakemake@output[['enrichment_figure_dir']])


########################################
## ANLYSIS
########################################

print(table(go_groups_dt[,go_group]))

## put some go term in the 'others' group
go_term_others <- unlist(as.data.table((table(go_groups_dt[,go_group])))[N <= 10, V1])
go_groups_dt[go_group %in% go_term_others, go_group := 'others']

## extrach the go groups into a vector
go_group_vector <- unique(unlist(go_groups_dt[,go_group]))


##########
## Generate barplot for each go group
##########


## initialize the list that contain the barchart
barchart_list <- list()

## extract the data that are not significant
go_groups_dt_without_ns <- go_groups_dt[!(data_type == 'ns'),]

## for each  go group, generate a barplot of the GO term with the pval
for (i_go_group_vector in seq(1, length(go_group_vector), 1)) {

    ## take the subset data that are associated with the go goup of interest
    subset <- copy(go_groups_dt_without_ns[go_group == go_group_vector[i_go_group_vector],])
    
    ## replace the space in the go group name by underscore
    go_group_name <- str_replace_all(
        go_group_vector[i_go_group_vector],
        pattern = ' ',
        replacement = '_'
    )

    ## initialize the file path
    file_path <- paste(
        snakemake@output[['enrichment_figure_dir']],
        '/barplot_',
        go_group_name,
        '.svg',
        sep = ''
    )

    ## reorder the go term based on the pvalue
    subset[, log_padjust := -log2(as.numeric(p.adjust))]
    subset <- subset[order(log_padjust, decreasing = F),]    
    #subset <- unique(subset[, c('data_type', 'Description', 'log_padjust'), with = F])
    subset$Description <- factor(
        subset$Description,
        levels = unique(unlist(subset$Description))
    )

    ## create the barplot
    GO_plot <- ggplot(
        data = subset,
        aes(
            x = Description,
            #y = -log2(as.numeric(p.adjust)),
            y = log_padjust,
            fill = data_type
        )
    ) +
        geom_bar(stat="identity") +
        coord_flip() +
        labs(
            title = go_group_vector[i_go_group_vector],
            y = '-Log10(pvalue)',
            x = 'Gene Ontology Term'
        ) +
        #scale_fill_manual(values = c('red3', 'blue3')) +
        scale_fill_manual(values = c('positive' = 'red3', 'negative' = 'blue3')) +
        theme_classic(base_size = 24)


    barchart_list[[go_group_name]] <- GO_plot

    ## save the figure into a file
    ggsave(
        plot = GO_plot,
        filename = file_path,
        device = 'svg',
        height = 24, 
        width = 20 #12
    )

}


print(barchart_list)


## arrange the plots together
barchart_all <- ggarrange(
    barchart_list[["establishment_of_localization"]],
    barchart_list[["cell_communication"]],
    barchart_list[["macromolecule_localization"]],
    barchart_list[["cellular_component_organization_or_biogenesis"]],
    barchart_list[["regulation_of_biological_quality"]],
    barchart_list[["system_process"]],
    nrow = 2,
    ncol = 3
    #heights = c(2,2,1)
    #widths = c(2,2,1)
)

svg(
    filename = paste(
        snakemake@output[['enrichment_figure_dir']],
        '/',
        "positive_correlation_goterm",
        '.svg',
        sep = ""
    ),
    height = 40,
    width = 50
)
print(barchart_all)
dev.off()


##########
## Diagramm de Sankey
##########


## calculate the ratio of go term associated with negative or positive correlation in a go group
go_groups_dt_without_ns <- go_groups_dt[!(data_type == 'ns'), ]

## initialize the data table that will contain the ratio information between negative and positive data
ratio_dt <- data.table(
    go_name = go_group_vector
)

for (i_go_group_vector in seq(1, length(go_group_vector), 1)) {
    subset <- go_groups_dt_without_ns[go_group == go_group_vector[i_go_group_vector],]

    ## counting the number of go term for positive and negative genes
    n_positive <- (nrow(subset[data_type == 'positive'])) + 0.0000001
    n_negative <- (nrow(subset[data_type == 'negative'])) + 0.0000001

    ## calculate the ratio between positive and negative    
    ratio_value <- n_positive/n_negative

    ## put the ratio value into the data table
    ratio_dt[go_name ==  go_group_vector[i_go_group_vector], ratio := ratio_value]

}

## order the data table with the ratio
ratio_dt <- ratio_dt[order(ratio),]

## order the go group with the ratio values
go_groups_dt$go_group <- factor(
    go_groups_dt$go_group,
    levels = unlist(ratio_dt[,go_name])
)

## generate the sankey diagram
correlation_sankey_diagram <- ggplot(
    go_groups_dt,
    aes(
        #y = Freq,
        axis1 = data_type,
        axis2 = go_group,
        label = data_type
    )
) +
    geom_alluvium(aes(fill = data_type)) +
    #geom_stratum(width = 1/12, fill = "grey", color = "black") +
    geom_stratum(alpha = 0, size = 2) +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 8) +
    scale_x_discrete(expand = c(.1, .1)) +
    scale_fill_manual(values = c('blue2', 'grey', 'red2')) +

    ggtitle("Gene Ontology enrichment of the correlated genes") +
    theme_classic(base_size = 24) +
    theme(legend.position = "none")


########################################
## Save plots
########################################

## save the sankey diagramm
ggsave(
    plot = correlation_sankey_diagram,
    filename = paste(
        snakemake@output[['enrichment_figure_dir']],
        '/sankey_diagram.svg',
        sep = ''
    ),
    device = 'svg',
    height = 15,
    width = 12
)


#Sys.sleep(10000)





#clinical_data <- fread(
#    '/run/media/hoangdong/linux/work/phd/projects/gaba_and_glutamate_pathways_in_glioma_v2/data/clinical_data/TCGA_clinical_data_cancer',
#    sep = ','
#)

#id <- fread(
#    '/run/media/hoangdong/linux/work/phd/projects/gaba_and_glutamate_pathways_in_glioma_v2/ensembl_id.csv',
#    sep = ','
#)

#expression_data <- fread(
#    '/run/media/hoangdong/linux/work/phd/projects/gaba_and_glutamate_pathways_in_glioma_v2/data/expression_data/TCGA_expression_data_fpkm_cancer',
#    sep = ','
#)[, genes := sapply(genes, function(x) str_split(x, pattern = '[.]')[[1]][1])]





#timer2 <- fread(
#    '/run/media/hoangdong/linux/work/phd/projects/gaba_and_glutamate_pathways_in_glioma_v2/data/timer2_data.csv',
#    sep = ','
#)

#print(colnames(clinical_data))


#sample_id_vector <- unlist(clinical_data[, Sample_ID])
#sample_id_vector <- sapply(sample_id_vector, function(x) str_replace_all(x, pattern = 'A$', replacement = ''))
#sample_id_vector <- sapply(sample_id_vector, function(x) str_replace_all(x, pattern = 'B$', replacement = ''))
#sample_id_vector <- sapply(sample_id_vector, function(x) str_replace_all(x, pattern = 'C$', replacement = ''))

##print(sample_id_vector)


##break
#timer2 <- timer2[cell_type %in% sample_id_vector,]



#print(dim(timer2))

##setdiff(sample_id_vector, unlist(timer2[, cell_type]))

#break



#gene_symbol <- map_id(
#    unlist(expression_data[,genes]),
#    id,
#    'gene_id',
#    'gene_name'
#)

#expression_data[, genes := gene_symbol]

#table(clinical_data[, Project_ID])


#sample_GBM <- clinical_data[Project_ID == 'TCGA-GBM', Case_ID]
#sample_LGG <- clinical_data[Project_ID == 'TCGA-LGG', Case_ID]
#expression_GBM <- expression_data[, c('genes', sample_GBM), with = F]
#expression_LGG <- expression_data[, c('genes', sample_LGG), with = F]




#create_fourchette <- function(
#    int_input
#) {
#    fourchette <- c(seq(1, int_input, 50), int_input)
#    return(fourchette)
#}


#test <- create_fourchette(ncol(expression_GBM))

#for (i in seq(1, (length(test) - 1), 1)) {

#    inferior_window <- test[i] + 1
#    superior_window <- test[i+1]

#    subset <- expression_GBM[, c(1, seq(inferior_window, superior_window, 1)), with = F]


#    file_name <- paste(
#        "GBM_",
#        inferior_window,
#        "_to_",
#        superior_window,
#        sep = ""
#    )

#    fwrite(
#        subset,
#        paste(
#            '/run/media/hoangdong/linux/work/phd/projects/gaba_and_glutamate_pathways_in_glioma_v2/sub_expression',
#            '/',
#            file_name,
#            sep = ''
#        ),
#        sep = ','
#    )

#    print(file_name)
#}


#break