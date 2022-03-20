

#######################################################
## GENERATE FIGURES from the epidemiological data
#######################################################

#######################################################
## Load the library
#######################################################

library(data.table)
library(tidyverse)
library(RColorBrewer)

#######################################################
## Load the data
#######################################################

## load the epidemiological data
epidemio_dt <- fread(
    snakemake@input[['subgroup_epidemiological_data_all']],
    sep = ','
)
epidemio_dt[cluster == 'IDHmut_CODEL', cluster := "IDHmut\nCODEL"]
epidemio_dt[cluster == 'IDHmut_nonCODEL', cluster := "IDHmut\nnonCODEL"]

epidemio_dt[epidemio_dt == ''] <- 'unknown'
epidemio_dt[is.na(epidemio_dt)] <- 'unknown'


## variable
variable_vector <- colnames(epidemio_dt)[!(colnames(epidemio_dt) %in% c('Case_ID', 'cluster', 'subcluster'))]

## qualitative variable names
qualitative_variables <- c(
    "gender",
    "Karnofsky Performance Score",
    "classification",
    "IDH status",
    "1p/19q codeletion",
    "Chr 7 gain/Chr 10 loss",
    "TERT promoter status",
    "MGMT promoter status",
    "ATRX status",
    "EGFR",
    "FGFR",
    "CDKN2A",
    "CDKN2B",
    "MYB",
    "MYBL1",
    "Grade",
    'GNB1','GABRD','ALDH4A1',
    "Histology"
)


#gaba_gluta_genes <- c(
#    'GNB1','GABRD','ALDH4A1','GRIK3','RIMKLA','SLC1A7','GNG12','PRKACB','GNG5','GNAI3','GLUL','PLA2G4A','GNG4','ADCY3','CAD','PPP3R1','GFPT1','KCNJ3','GAD1','GLS','PLCL1','TRAK2','CPS1','ITPR1','GRM7','SLC6A11','SLC6A1','SLC38A3','GNAI2','GRM2','CACNA1D','NIT2','ADCY5','TRPC1','PLD1','GNB4','NAT8L','GABRG1','GABRA2','GABRA4','GABRB1','PPAT','PPP3CA','GRIA2','ADCY2','SLC1A3','HOMER1','GRIA1','GABRB2','GABRA1','GABRG2','GFPT2','ALDH5A1','GABBR1','ITPR3','GRM4','GRIK2','DDO','GRM1','ADCY1','ASL','GNAI1','GRM3','GNG11','ASNS','GNB2','GRM8','PPP3CC','ADCY8','GPT','SLC1A1','GNAQ','GABBR2','GRIN3A','GNG10','ASS1','GRIN1','CACNA1B','GAD2','PPP3CB','GLUD1','GOT1','SLC17A6','SLC1A2','FOLH1','ASRGL1','GNG3','PLCB3','GRK2','SHANK2','GRM5','GRIA4','GRIK4','SLC6A13','SLC6A12','CACNA1C','GNB3','RIMKLB','GABARAPL1','GRIN2B','ITPR2','SLC38A1','SLC38A2','ADCY6','SLC17A8','ADCY4','GNG2','GPHN','GABRB3','GABRA5','GABRG3','PLCB2','JMJD7-PLA2G4B','GNB5','HOMER2','ADCY9','ABAT','GRIN2A','PRKCB','MAPK3','GPT2','ADCY7','GNAO1','GOT2','GABARAPL2','ASPA','PLD2','DLG4','GABARAP','HAP1','NSF','PRKCA','GRIN2C','DLGAP1','GNG7','CACNA1A','PRKACA','SLC1A6','HOMER3','GRIK5','PLA2G4C','GRIN2D','SLC17A7','IL4I1','SHANK1','PRKCG','PLCB1','PLCB4','SRC','SLC32A1','SLC12A5','GNAS','GRIK1','KCNJ6','MAPK1','GRK3','ADSL','SLC38A5','GLUD2','GRIA3','GABRE','GABRA3','GABRQ','MYB','EGFR','FGFR','MYBL1','CDKN2A','CDKN2B'
#)

#qualitative_variables <- c(
#    qualitative_variables,
#    gaba_gluta_genes
#)


## qualitative variable names
quantitative_variables <- c(
    "age_at_index",
    "TumorPurity"
)

## vector of the cluster names
cluster_vector <- (unique(unlist(epidemio_dt[, cluster])))

## vector of the subcluster names
subcluster_vector <- unique(unlist(epidemio_dt[, classification2]))




## colors that will be used for the generation of the plots
color_palette <- brewer.pal(9, "Set1")

print(color_palette[1:6])

## create a datatable containing the association of colors for each cluster name
color_data <- data.table(
    cluster = sort(cluster_vector),
    color = color_palette[1:(length(cluster_vector))]
)

## put the orange color for the mixed cluster
color_data[cluster == "MIXED", color :=  "#FF7F00"]



#######################################################
## Fonctions for figure generation
#######################################################

## function for creating the plot
generate_qualitative_plot <- function(
    data_input,
    variable_name,
    title_name_input,
    file_path
) {

    ## remove the underscore in the title name
    title_name_input <- str_replace_all(
        string = title_name_input,
        pattern = "_",
        replacement = " "
    )

    ## generate the color palette
    color_palette <- c(
        brewer.pal(9, "Set1"),
        brewer.pal(8, "Dark2")
    )

    ## generate the plot
    plot <- ggplot(
        data = data_input,
        aes(
            x = cluster,
            y = percentage,
            fill = variable
        )
    ) +
        geom_bar(
            position = "fill",
            stat = "identity"
        ) +
        labs(
            title = title_name_input,
            x = "Cluster",
            y = "Percentage",
            fill = variable_name
        ) +
        theme_classic(
            base_size = 32
        ) +
        scale_fill_manual(
            values = color_palette # color
        ) 

    ## save the plot
    plot <- ggsave(
        plot = plot,
        filename = file_path,
        device = "svg",
        width = 16,
        height = 10
    )

    ## return the plot
    return(plot)
}

## function for creating the plot using the quantitative data
generate_quantitative_plot <- function(
    data_input,
    variable_name,
    title_name_input,
    file_path
) {

    ## generate the color palette
    color_palette <- c(
        brewer.pal(9, "Set1"),
        brewer.pal(8, "Dark2")
    )
 

    ## generate the plot
    plot <- ggplot(
        data = data_input,
        aes(
            x = as.factor(cluster),
            y = as.numeric(unlist(data_input[, variable_name, with = F])),
            fill = as.factor(cluster)
        )
    ) +
        geom_boxplot(
            outlier.colour="black",
            outlier.shape=16,
            outlier.size=2,
            notch = TRUE,
            lwd = 4
        ) +  
        labs(
            title = title_name_input,
            x = "Cluster",
            y = variable_name,
            fill = variable_name
        ) +
        theme_classic(
            base_size = 32
        ) +
        theme(legend.position = "none") +
        scale_fill_manual(
            values = unlist(color_data[,color]) # color
        ) 

    ## save the plot
    ggsave(
        plot = plot,
        filename = file_path,
        device = "svg",
        width = 16,
        height = 10
    )

    ## return the plot
    return(plot)
}


#######################################################
## Generate figures for each cluster
#######################################################

## create the directory that will contain the figures
dir.create(snakemake@output[['generated_figures']])

for (i_variable in seq(1, length(variable_vector), 1)) {
    
    ## variable name
    variable_name <- variable_vector[i_variable]

    ## formate the variable name to put it in a path file
    variable_name2 <- str_replace_all(
        variable_name,
        pattern = '[/]|[ ]',
        replacement = '_'
    )

    ## path of the figure file
    figure_path <- paste(
        snakemake@output[['generated_figures']],
        '/',
        variable_name2,
        '.svg',
        sep = ''
    )
    #print(figure_path)

    ## extract the data associated with the variable name
    subdata <- copy(epidemio_dt[,c('cluster', variable_name), with = F])

    ## look the number of categories of the variable for knowing if the variable is quantitative of qualitative
    nb_categories <- dim(table(subdata[, variable_name, with = F]))
    
    if (variable_name %in% qualitative_variables) { ## we take the variable as a category variable

        ## add some information such as the count and percentage
        subdata2 <- unique(subdata[, count := (.N), by = c('cluster', variable_name)])
        subdata2[, total := sum(count), by = c('cluster')]
        subdata2[, percentage := round((count/total*100), digit = 2), by = c('cluster')]
        subdata2 <- subdata2[, setnames(.SD, variable_name, 'variable')] # rename the column

        if (variable_name == 'Karnofsky Performance Score') {

            subdata2$variable <- factor(
                subdata2$variable,
                levels = c('100', '90','80','70','60','50','40','unknown')
            )
        }


        # generate the figure
        plot <- generate_qualitative_plot(
            data_input = subdata2,
            variable_name = variable_name,
            title_name_input = variable_name,
            file_path = figure_path
        )

    }
    else { ## else the variable is quantitative
        
        plot <- generate_quantitative_plot(
            data_input = subdata,
            variable_name = variable_name,
            title_name_input = variable_name,
            file_path = figure_path
        )
    }
}

#Sys.sleep(10000)



########################################################
### Generate figures for each subcluster
########################################################


## create the directory that will contain the figures of the subcluster
subcluster_dir <- paste(
    snakemake@output[['generated_figures']],
    "/subcluster",
    sep = ''
)

dir.create(
    subcluster_dir,
    recursive = T
)

for (i_variable in seq(1, length(variable_vector), 1)) {

    for (j_subcluster in seq(1, length(subcluster_vector), 1)) {
        
        ## extract the variable, cluster and subcluster name of interest
        variable_name <- variable_vector[i_variable]
        #cluster_name <- cluster_vector[j_cluster]
        subcluster_name <- subcluster_vector[j_subcluster]

        ## formate the variable name to put it in a path file
        variable_name2 <- str_replace_all(
            variable_name,
            pattern = '[/]|[ ]',
            replacement = '_'
        )

        ## path of the figure file
        figure_path <- paste(
            subcluster_dir,
            '/',
            variable_name2,
            '_',
            subcluster_name,
            '.svg',
            sep = ''
        )

        print(figure_path)

        ## extract the data associated with the variable name and the subcluster name
        subdata <- copy(epidemio_dt[
            classification2 == subcluster_name,
            c('cluster', variable_name),
            with = F
        ])

        if (variable_name %in% qualitative_variables) { 
            ## add some information such as the count and percentage
            subdata2 <- unique(subdata[, count := (.N), by = c('cluster', variable_name)])
            subdata2[, total := sum(count), by = c('cluster')]
            subdata2[, percentage := round((count/total*100), digit = 2), by = c('cluster')]
            subdata2 <- subdata2[, setnames(.SD, variable_name, 'variable')] # rename the column

            ## generate the figure
            plot <- generate_qualitative_plot(
                data_input = subdata2,
                variable_name = variable_name,
                title_name_input = subcluster_name,
                file_path = figure_path
            )

        }
        else { ## else the variable is quantitative

            ## generate the figure of the quantitative plot
            plot <- generate_quantitative_plot(
                data_input = subdata,
                variable_name = variable_name,
                title_name_input = subcluster_name,
                file_path = figure_path
            )
        }
    }
}
