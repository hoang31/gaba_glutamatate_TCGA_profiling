
###################################################
###################################################


##### The scripts contain all the general functions


###################################################
###################################################


##### load library
library(data.table)


###################################################
###################################################


###### Function for transform the ensembl ID to gene symbols

##### inputs

## vector of ensembl id that you want to transform to gene symbols
## id table with two columns : "gene_id" (corresponding to ensembl id) and "gene_name" (corresponding to the gene symbols)

ensemblID_to_geneSymbole <- function(ensemblID, id_table) {

    ## remove the point
    ensemblID <- sapply(ensemblID, function(x) strsplit(x, split = "[.]")[[1]][1])

    ## extract the ids which does not exist in the id table
    differences <- setdiff(ensemblID, unlist(id_table[, gene_id]))
    #differences <- setdiff(unlist(id_table[, gene_id]), ensemblID)

    dt_diffrences <- as.data.table(cbind(differences, differences)) # transform to data table

    ## rename the colnames
    colnames(dt_diffrences) <- colnames(id_table)

    ## add the difference data table into the id table
    id_table <- rbind(id_table, dt_diffrences)
    
    ### set indexes in the data table
    #setkey(id_table, gene_id)

    ## extract the gene name
    output <- id_table[gene_id %in% ensemblID, ]


    ## reorder with the ensemblID
    output_reordered <- output[order(match(output[, gene_id], ensemblID))]

    return(output_reordered[, gene_name])
}


###################################################
###################################################


##### function for retreinve the directory path and add a name file into it

##### inputs

## PATH_INPUT : path of the file ; example : "/home/xxx/xxx"
## FILE_NAME_INPUT : file name ; example : "file_name.txt"

create_file_path <- function(PATH_INPUT, FILE_NAME_INPUT) {
    file_name_path_output <- paste(PATH_INPUT,
        FILE_NAME_INPUT,
        sep = "/"
    )
    return(file_name_path_output)
}


###################################################
###################################################


## create the survival curve from a object fitted by the R packages survminer and survival
create_survival_curve <- function(
    FIT_MODEL_INPUT, # survival data containing 4 columns : sample_id, cluster, status (0 = alive, 1 = dead), time (in months),
    FILE_PATH_INPUT, # file path of the generated survival curve
    COLOR_PALETTE_INPUT,
    DISPLAY_PVALUE = FALSE, # dosplay the pvalue
    WIDTH_INPUT,
    HEIGTH_INPUT

) {
    ## if there is a path for saving the file
    if (!(missing(FILE_PATH_INPUT))) {
        generate_plot <- TRUE
    }

    if (!(missing(FILE_PATH_INPUT))) {
        generate_plot <- TRUE
    }
    
    if ((missing(COLOR_PALETTE_INPUT))) {
        ## choose color
        color_palette <- brewer.pal(9, "Set1")

        ## number of color needed for the survival plot
        n_col <- length(unique(unlist(survival_data[, "cluster", with = F])))

        ## extract the colors
        COLOR_PALETTE_INPUT <- color_palette[1:n_col]
    }
    
    ## fit the survival data
    # fit <- survfit(Surv(time, status) ~ cluster, data = SURVIVAL_DATA_INPUT)

    ## create the survival plot
    survival_plot <- ggsurvplot(
        FIT_MODEL_INPUT,
        pval = DISPLAY_PVALUE, 
        # pval = TRUE,
        surv.median.line = "hv",
        conf.int = FALSE,
        #risk.table = TRUE,
        #tables.height = 0.2,
        #tables.theme = theme_cleantable(),
        # break.time.by = 10,
        # xlim = c(0,60),
        #legend.title = "Glioma Cluster",
        legend.title = "",
        palette = COLOR_PALETTE_INPUT,
        size = 4

    ) + 
    labs(
        title = "",
        x = "Time (Months)"
    )

    survival_plot$plot <- survival_plot$plot + 
        theme_bw(base_size = 40) +
        theme(legend.position="top")


    ## save the survival plot
    svg(
        filename = FILE_PATH_INPUT,
        width = WIDTH_INPUT,
        height = HEIGTH_INPUT
    )
    print(survival_plot)
    dev.off()
}




###################################
## function for change a id type to a other
###################################

map_id <- function(
    id_vector, # vector of ids that have to be change
    id_table, # table that contain two column of id types
    from_id_type, # id type that will be read (correspond to the name of one of the column of the id_table)
    to_id_type # id type that will be changed to (correspond to the name of the other column of the id_table)
) {

    ## extract from the id_table all the information related to the id_vector
    id_of_interest_dt <- id_table[
        unlist(id_table[,from_id_type, with = F]) %in% id_vector,
    ]

    ## extract the id that are not in the id_table
    if (length(id_vector) > nrow(id_table)) {
        lacked_id <- setdiff(id_vector, unlist(id_table[, from_id_type, with = F]))

        lacked_dt <- data.table(
            v1 = lacked_id,
            v2 = lacked_id
        )

        colnames(lacked_dt) <- colnames(id_table)
    }

    ## extract the id that are not in the id_table
    if (length(id_vector) < nrow(id_table)) {
        lacked_id <- setdiff(id_vector, unlist(id_table[, from_id_type, with = F]))

        lacked_dt <- data.table(
            v1 = lacked_id,
            v2 = lacked_id
        )

        colnames(lacked_dt) <- colnames(id_table)
    }

    ## rbind with the data that are lacking
    id_of_interest_dt <- rbind(
        id_of_interest_dt,
        lacked_dt
    )

    ## reorder the data table depending of the order of the id vector
    id_of_interest_dt <- id_of_interest_dt[order(match(unlist(id_of_interest_dt[, from_id_type, with = F]), id_vector))]

    ## extrract the id type of interest
    id_of_interest <- unlist(id_of_interest_dt[, to_id_type, with = F])

    ## remove names
    names(id_of_interest) <- NULL

    ## return the output
    return(id_of_interest)

}



###################################
## function for performing the data table transposition keeping names
###################################

transpose_datatable <- function(
    data_table, # data frame or data table
    column_name = NULL, # column name that correspond to the final column name
    new_name = "V1" # character that correspond to the column name that will contain the values of the old column names
) {

    ## copy the data table
    transposed_dt <- as.matrix(copy(data_table))

    ## do the transposition
    transposed_dt <- as.data.table(t(transposed_dt), keep.rownames = T)

    ## delete the row that correspond to the colnames
    if (is_null(column_name) == FALSE) {

        ## extract the column names into a vector
        column_name_vector <- unlist(transposed_dt[rn == column_name, ])
        
        ## delete the row that correspond to the column names
        transposed_dt <- transposed_dt[rn != column_name, ]

        ## rename column 
        colnames(transposed_dt) <- column_name_vector

        ## rename the column related to the previous column names
        transposed_dt <- transposed_dt[, setnames(.SD, column_name, new_name)] 
    }

    return(transposed_dt)

}



###################################
## function for performing pair statistical test
###################################

## Function for performing the statistical analysis between each groups
statistical_analysis_qualitative <- function(
    data_table_input, # data table that contain in row the samples and in column the different variables
    variable_name, # variable name that correspond to the column of interest
    variable_grouping_name # variable name that correspond to the grouping variable
) {
    
    ## extract the data associated with the variable name and the cluster columns
    data_table_input <- copy(data_table_input[, c(variable_name, variable_grouping_name), with = F])

    ## generate the contingency table associated with the variable of interest and the cluster
    contengency_table <- as.data.table(as.data.frame.matrix(table(data_table_input)), keep.rownames = "category")

    ## initialize the data table that will contain the statistical results
    statistical_results_dt <- data.table(
        "variable" = character(),
        "cluster1" = character(),
        "cluster2" = character(),
        "pvalue" = numeric()
    )

    ## identification of the cluster names
    cluster_name_vector <- unique(unlist(data_table_input[, cluster]))

    ## do the paired combination of the cluster name vector
    cluster_combination <- combn(
        x = cluster_name_vector,
        m = 2
    )
        
    for (i_combination in seq(1, ncol(cluster_combination), 1)) {
        
        ## extract the clusters of interest
        cluster1 <- cluster_combination[1, i_combination]
        cluster2 <- cluster_combination[2, i_combination]

        #print("=================")
        
        ## extract the contingency data associated with
        subset_contingency_table <- copy(contengency_table[, c(cluster1, cluster2), with = F ])
        #print(subset_contingency_table)

        ## do the statistical test
        pval <- fisher.test(
            x = subset_contingency_table,
            simulate.p.value = TRUE
        )$p.value

        ## put the results into a data table
        results_dt <- data.table(
            "variable" = variable_name,
            "cluster1" = cluster1,
            "cluster2" = cluster2,
            "pvalue" = pval
        )

        ## put the results into the statistical_results_dt
        statistical_results_dt <- rbind(
            statistical_results_dt,
            results_dt
        )

    }
        
    

    ## return the statistical results
    return(statistical_results_dt)

}


###################################
## function for performing pair statistical test for quantitative variable
###################################

statistical_analysis_quantitative <- function(
    data_table_input, # data table that contain in row the samples and in column the different variables
    variable_name, # variable name that correspond to the column of interest
    variable_grouping_name # variable name that correspond to the grouping variable
) {
    
    ## extract the data associated with the variable name and the cluster columns
    data_table_input <- copy(data_table_input[, c(variable_grouping_name, variable_name), with = F])

    ## initialize the data table that will contain the statistical results
    statistical_results_dt <- data.table(
        "variable" = character(),
        "cluster1" = character(),
        "cluster2" = character(),
        "pvalue" = numeric()
    )

    ## identification of the cluster names
    cluster_name_vector <- unique(unlist(data_table_input[, cluster]))

    ## do the paired combination of the cluster name vector
    cluster_combination <- combn(
        x = cluster_name_vector,
        m = 2
    )
        
    for (i_combination in seq(1, ncol(cluster_combination), 1)) {
        
        ## extract the clusters of interest
        cluster1 <- cluster_combination[1, i_combination]
        cluster2 <- cluster_combination[2, i_combination]

        ## retreive the values related to each cluster and cell type of interest
        cluster1_values <- as.numeric(unlist(data_table_input[cluster == cluster1, variable_name, with = F]))
        cluster2_values <- as.numeric(unlist(data_table_input[cluster == cluster2, variable_name, with = F]))

        ## do the statistical test
        pval <- wilcox.test(
            cluster1_values,
            cluster2_values
        )$p.value

        ## put the results into a data table
        results_dt <- data.table(
            "variable" = variable_name,
            "cluster1" = cluster1,
            "cluster2" = cluster2,
            "pvalue" = pval
        )

        ## put the results into the statistical_results_dt
        statistical_results_dt <- rbind(
            statistical_results_dt,
            results_dt
        )

    }
        
    

    ## return the statistical results
    return(statistical_results_dt)

}


###################################
## function for getting the percentage of a categorical vcariable
###################################

get_percentage_qualitative <- function(
    data_table_input, # data table that contains in column the clinical variables and in rowns each samples
    variable_name_input, # name of the column that corresponds to the variable of interest
    group_name_input, # name of the column that corresponds the grouping variable
    formated = F # formate the output without melting the data
) {

    ## get the group name into a vector
    group_name_vector <- sort(unique(unlist(data_table_input[, group_name_input, with = F])))

    ## retreive the data associated with the variable and the group
    data_table_input <- data_table_input[, c(variable_name_input, group_name_input), with = F]

    ## retreive the number of sample
    sample_nb <- nrow(data_table_input)

    ## do the counting
    counting_table <- as.data.frame.matrix(table(data_table_input))

    ## transform the data frame to data table
    setDT(counting_table, keep.rownames = "category")

    ## copy the counting table
    counting_table_percentage <- copy(counting_table)

    ## do the percentage
    counting_table_percentage <- counting_table_percentage[
        ,
        (group_name_vector) := lapply(
            .SD,
            function(x) round(x/sum(x)*100, digit = 2)
        ),
        .SDcols = group_name_vector
    ]
    ## copy the counting table
    counting_table_percentage_formated <- copy(counting_table)

    ## do the percentage reformating the data
    counting_table_percentage_formated <- counting_table_percentage_formated[
        ,
        (group_name_vector) := lapply(
            .SD,
            function(x) {
                
                ## calculate the percentage
                percentage <- round(x/sum(x)*100, digit = 2)

                ## generate the caracters
                output <- paste(
                    x,
                    "/",
                    sum(x),
                    " (",
                    percentage,
                    "%)",
                    sep = ""
                )

                return(output)
            }
        ),
        .SDcols = group_name_vector
    ]

    ## put the variable name into the counting table
    counting_table_percentage_formated[, variable := variable_name_input]

    ## if formated == TRUE, we will return the countint table before melting the data    
    if (formated == TRUE) {
        return(counting_table_percentage_formated)
    }

    ## melt the data
    counting_table_percentage_melted <- melt(
        counting_table_percentage,
        id = 'category', 
        measure = group_name_vector
    )

    ## put factor for the group variable
    counting_table_percentage_melted$variable <- factor(
        counting_table_percentage_melted$variable,
        levels = group_name_vector
    )

    ## return the counting table percentage
    return(counting_table_percentage_melted)
}


###################################
## function for getting the percentage of a categorical vcariable
###################################

get_percentage_quantitative <- function(
    data_table_input, # data table that contains in column the clinical variables and in rowns each samples
    variable_name_input, # name of the column that corresponds to the variable of interest
    group_name_input, # name of the column that corresponds the grouping variable
    formated = F # formate the output without melting the data
) {

    ## get the group name into a vector
    group_name_vector <- sort(unique(unlist(data_table_input[, group_name_input, with = F])))

    ## retreive the data associated with the variable and the group
    data_table_input <- data_table_input[, c(variable_name_input, group_name_input), with = F]

    ## retreive the number of sample
    sample_nb <- nrow(data_table_input)

    ## get all the group value into a vector
    group_name_vector <- sort(unique(unlist(data_table_input[, group_name_input, with = F])))

    ## for each cluster, retreive the summary of the value
    summary_by_cluster_list <- lapply(
        X = group_name_vector,
        function(x) {

            ## get the position of samples associated with the group name of interest
            pos <- grep(x = unlist(data_table_input[,group_name_input, with = F]), pattern = x)

            ## extract the value associated with the group name samples
            summary_values <- as.data.table(as.matrix(summary(as.numeric(unlist(data_table_input[pos, variable_name_input, with = F])))), keep.rownames = "category")[, setnames(.SD, 'V1', x)]

            return(summary_values)
        }
    )

    ## merge all the sublists into a data table
    summary_by_cluster_dt <- Reduce(
        x = summary_by_cluster_list,
        function(x, y) {
            merge(
                x,
                y,
                by = "category",
                sort = F
            )
        }
    )

    ## add the variable into the summary cluster data table
    summary_by_cluster_dt[, variable := variable_name_input]

    ## if formated == TRUE, we will return the countint table before melting the data    
    if (formated == TRUE) {
        return(summary_by_cluster_dt)
    }
    else {
        return(data_table_input)
    }

}
