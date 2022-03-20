##################################################
###################################################


## Differential gene expression analysis


###################################################
###################################################


## LOAD LIBRARIES
library(data.table)
library(stringr)


###################################################
###################################################


##### LOAD DATA

## load clinical data
d_clinical <- fread(
    file = snakemake@input[["clinical_data"]],
    sep = ",",
    header = T
)

## load the external TCGA clinical data
d_clinical_external <- fread(
    file = snakemake@input[["external_clinical_data"]],
    sep = ","
) 

## load the sample clusters
d_clusters <- fread(
    file = snakemake@input[["cluster_group"]],
    sep = ","
)

## load the subcluster data
#d_subclusters <- fread(
#    file = snakemake@input[["subcluster_group"]],
#    sep = ","
#)

## load the purity data
d_purity <- fread(
    file = snakemake@input[["purity_score_tumor"]],
    sep = ","
)

## load the alteration data
d_alteration <- fread(
    file = snakemake@input[["cnv_data"]],
    sep = ","
)

cat("\n\n##################################################")
cluster_counting <- table(d_clusters[,"cluster"])
print(cluster_counting)
cat("###################################################\n\n")

## load external scripts
source(snakemake@params[["utils_R"]])











###################################################
###################################################

# d_clinical <- fread(
#     file = "~/Documents/gaba_and_glutamate_pathways_in_glioma/data/clinical_data/TCGA_clinical_data_cancer",
#     sep = ",",
#     header = T
# )

# ## load the external TCGA clinical data
# d_clinical_external <- fread(
#     file = "~/Documents/gaba_and_glutamate_pathways_in_glioma/data/clinical_data/TCGA_clinical_EXTERNAL.csv",
#     sep = ","
# ) 


# ## load the external TCGA clinical data
# d_TCGA_cluster <- fread(
#     file = "~/Documents/TCGA_clusters.csv",
#     sep = ","
# ) 




# d_clinical2 <- copy(d_clinical)



# ## extract the id for the count files
# d_clinical2 <- unique(
#     d_clinical2[
#         grep(
#             x = d_clinical2[, File_Name],
#             pattern = "counts"
#         ),
#     ]
# )


# names(d_clinical2)

# ## remove duplicate column (Project ID)
# d_clinical2 <- d_clinical2[
#     ,
#     unique(names(d_clinical2)),
#     with = FALSE
# ]
# ## merge clinical data with other clinical data
# d_clinical2 <- merge(
#     d_clinical2,
#     d_clinical_external,
#     by.x = "Case_ID",
#     by.y = "Case",
#     all.x = T,
#     sort = F
# )

# d_clinical_IDH_wt <- d_clinical2[`IDH status` == "WT", "Case_ID"]

# ## NORMAL LIKE glioma previously identified
# NL_samples <- c(
#     "TCGA-DU-A76K",
#     "TCGA-P5-A5EY",
#     "TCGA-DB-A75P",
#     "TCGA-HT-8558",
#     "TCGA-CS-6669",
#     "TCGA-DU-7292",
#     "TCGA-QH-A6XC",
#     "TCGA-HT-8015",
#     "TCGA-HT-8564",
#     "TCGA-FG-8181",
#     "TCGA-HT-8107",
#     "TCGA-HT-8019",
#     "TCGA-DU-8162",
#     "TCGA-FG-7643"
# )


# d_clinical_IDH_wt[Case_ID %in% NL_samples, "cancer type" := "NL" ]
# d_clinical_IDH_wt[!(Case_ID %in% NL_samples), "cancer type" := "OT" ]


# d_clusters <- d_clinical_IDH_wt[, c("Case_ID", "cancer type")][, setnames(.SD, c("Case_ID", "cancer type"), c("sample_id", "cluster"))]
# cluster_counting <- table(d_clusters[,"cluster"])
# print(colnames(d_clinical_external))

# fwrite(
#     d_clusters[order(cluster),],
#     file = "~/Documents/normal_like_glioma_project/project/data/d_clusters.csv",
#     sep = ","
# )
# break
###################################################
###################################################




















###################################################
###################################################


##### FUNCTION FOR EXTRACTING SPECIFIC DATA OF EACH CLUSTER
extract_data <- function(
    CLUSTER_GROUP_INPUT, # data table of two columns ( "sample_id" and "cluster")
    CLINICAL_DATA_INPUT, # clinical data table : columns corresponding to clinical (+ one column called "Gen_Name corresponding to the gene names"), rows corresponding to genes,
    QUALITATIVE_VARIABLE_INPUT,
    QUANTITATIVE_VARIABLE_INPUT

) {
    
    ## extract the variables of interest
    CLINICAL_DATA_INPUT <- CLINICAL_DATA_INPUT[, c( "Case_ID",QUALITATIVE_VARIABLE_INPUT, QUANTITATIVE_VARIABLE_INPUT), with = F]

    ## replace the na value to unknown
    CLINICAL_DATA_INPUT[is.na(CLINICAL_DATA_INPUT)] <- "unknown"

    ## initialize the epidemiological data associated with each clusters
    cluster_epidemiological_data_list <- list()

    ## extract the cluster names from the cluster group data
    cluster_names <- unique(unlist(CLUSTER_GROUP_INPUT[, "cluster"]))

    ## for each cluster, extract the N first most expressed
    for (i_cluster_name in seq(1, length(cluster_names), 1)) {
        
        # print("----------------------------")
        # print(cluster_names[i_cluster_name])

        ## extract the sample ids associated with the cluster group
        cluster_sample <- unlist(CLUSTER_GROUP_INPUT[cluster == cluster_names[i_cluster_name], "sample_id", with = F])

        ## extract the clinical data associated wutg the cluster sample
        cluster_clinical_data <- CLINICAL_DATA_INPUT[Case_ID %in% cluster_sample,]

        # print("----------------------")
        # print(cluster_names[i_cluster_name])
        # print("----------------------")

        # print(names(cluster_epidemiological_data_list))
        for (i_variable_name in seq(2, ncol(cluster_clinical_data))) {
            # print(cluster_epidemiological_data_list)

            ## retreive the variable name
            variable_name <- colnames(cluster_clinical_data)[i_variable_name]

            if (variable_name %in% QUALITATIVE_VARIABLE_INPUT) {

                ## counting table
                counting_table <- as.table(table(cluster_clinical_data[, variable_name, with = F]))

                # print(counting_table)
                # cat("\n##########################\n")
                # print(variable_name)
                # print(cluster_names[i_cluster_name])
                # print(counting_table)
                # print(sum(counting_table))
                # cat("\n##########################\n")
                # print(cluster_epidemiological_data_list)
                # print(cluster_epidemiological_data_list)

                ## add the information into the list
                cluster_epidemiological_data_list[[variable_name]][[cluster_names[i_cluster_name]]] <- counting_table

            }

            if (variable_name %in% QUANTITATIVE_VARIABLE_INPUT) {
                
                ## extract the age as a numeric vector
                values_vector <- as.numeric(unlist(cluster_clinical_data[, variable_name, with = F]))

                ## remove the NA values
                values_vector <- values_vector[!is.na(values_vector)]

                cluster_epidemiological_data_list[[variable_name]][[cluster_names[i_cluster_name]]] <- values_vector

            }
        }
    }

    ## return the output
    return(cluster_epidemiological_data_list)


}

###################################################
###################################################


##### FUNCTION FOR PERFORMING THE STATISTICAL TEST FOR EACH VARIABLES
statistic <- function(
    EPIDEMIO_DATA_LIST_INPUT,
    QUALITATIVE_VARIABLE_INPUT,
    QUANTITATIVE_VARIABLE_INPUT
) {

    ## extract cluster names
    variable_names <- names(EPIDEMIO_DATA_LIST_INPUT)
    # variable_names <- variable_names[variable_names != "gender"]

    ## extract variables
    cluster_names <- names(EPIDEMIO_DATA_LIST_INPUT[[1]])

    ## intialize the list of data 
    output_list <- list()

    
    for (i in variable_names) {
        
        print(i)

        ## data table that will contains the pvalue with the clusters used for the statistical test
        dt_pvalues <- data.table(
            "cluster1" = character(),
            "cluster2" = character(),
            "pvalue" = numeric()
        )

        ## if the variable are qunantitative, to a wilcoxon test
        if (i %in% QUANTITATIVE_VARIABLE_INPUT) {
            # cat("\n################################\n")
            # print(i)
            # print("quantitative")

            ## initialize the summary table
            summary_table <- data.table(
                "summary" = c(
                    "Min.",
                    "1st Qu.",
                    "Median",
                    "Mean",
                    "3rd Qu.", 
                    "Max."
                )
            )

            ## extract the data associated with the variable
            variable_data <- EPIDEMIO_DATA_LIST_INPUT[[i]]

            ## do all the combination for the statistical test
            combinations <- combn(
                x = seq(1, length(cluster_names), 1),
                m = 2
            )

            for (combination_i in seq(1, ncol(combinations), 1)) {

                # cat("\n-------\n")

                ## extract the cluster name and the values associated with the cluster
                cluster_name1 <- (names(variable_data[combinations[1,combination_i]]))
                v1 <- unlist((variable_data[combinations[1,combination_i]]))

                ## extract the cluster name and the values associated with the cluster
                cluster_name2 <- (names(variable_data[combinations[2,combination_i]]))
                v2 <- unlist((variable_data[combinations[2,combination_i]]))

                ## merge the summary table with the summary table associated with the cluster name 1
                if (!(cluster_name1 %in% colnames(summary_table))) {

                    ## retreive the summary of the cluster name1 values and rename the columns
                    cluster_summary_table <- as.data.table((as.matrix(summary(v1))), keep.rownames = T)[, setnames(.SD, c("rn", "V1"), c("summary", cluster_name1))]

                    ## merge the this cluster summary table with the main summary table
                    summary_table <- merge(
                        x = summary_table,
                        y = cluster_summary_table,
                        by = "summary",
                        all.x = T,
                        sort = F
                    )
                }

                ## merge the summary table with the summary table associated with the cluster name 2
                if (!(cluster_name2 %in% colnames(summary_table))) {

                    ## retreive the summary of the cluster name1 values and rename the columns
                    cluster_summary_table <- as.data.table((as.matrix(summary(v2))), keep.rownames = T)[, setnames(.SD, c("rn", "V1"), c("summary", cluster_name2))]

                    ## merge the this cluster summary table with the main summary table
                    summary_table <- merge(
                        x = summary_table,
                        y = cluster_summary_table,
                        by = "summary",
                        all.x = T,
                        sort = F
                    )
                }

                ## do the statistical test
                statitiscal_results <- wilcox.test(
                    x = v1,
                    y = v2
                )

                ## extract the pvalue and the clusters names used for the stat test
                pvalue <- as.data.table(
                    list(
                        cluster_name1,
                        cluster_name2,
                        statitiscal_results[["p.value"]]
                    ),
                )[, setnames(
                    .SD,
                    c("V1", "V2", "V3"),
                    colnames(dt_pvalues)
                )]
                
                ## add the pvalue to the pvalues data table 
                dt_pvalues <- rbind(
                    dt_pvalues,
                    pvalue
                )

            }

            ## put the summary table into the output list
            output_list[[i]][["data"]] <- summary_table
            output_list[[i]][["p_values"]] <- dt_pvalues

            # print("ca bug pas1")
        }
        

        ## if the variable is qualitative, do a fisher test/chi2 test
        if (i %in% QUALITATIVE_VARIABLE_INPUT) {

            cat("\n################################\n")
            # print(i)
            # print("qualitative")
            # print("ca bug pas2")
            ## extract the data associated with the variable
            variable_data <- EPIDEMIO_DATA_LIST_INPUT[[i]]

            var <- c()
            for (i_cluster_name in seq(1, length(variable_data), 1)) {

                ## retreive the data 
                cluster_data <- variable_data[[cluster_names[i_cluster_name]]]
                
                ## put "unknown"
                if (length(names(cluster_data)) >= 3 ) {
                    names(variable_data[[cluster_names[i_cluster_name]]])[names(variable_data[[cluster_names[i_cluster_name]]]) == ""] <- "unknown"
                    names(cluster_data)[names(cluster_data) == ""] <- "unknown"

                }

                ## save all the variables
                var <- append(var, names(cluster_data))
            }

            var <- unique(var)

            ## initialize the data that will contain the data table with the variable
            dt_data <- data.table(
                "variable" = var

            )
            ## add data associated with each cluster into the dt_data
            for (i_cluster_name in seq(1, length(variable_data), 1)) {

                ## retreive the cluster name
                cluster_name1 <- cluster_names[i_cluster_name]

                ## retreive the data, transform it to data table and change the column names
                cluster_data <- as.data.table(
                    variable_data[[cluster_names[i_cluster_name]]]
                )[
                    ,
                    setnames(
                        .SD,
                        c("V1", "N"),
                        c("variable", cluster_name1)
                    )
                ]

                ## merge the cluster data with the dt_data
                dt_data <- merge(
                    dt_data,
                    cluster_data,
                    by = "variable",
                    all.x = T,
                    sort = F
                )

            }


            ## replace the NA values to 0
            dt_data[is.na(dt_data)] <- 0

            ## do all the combination for the statistical test
            combinations <- combn(
                x = seq(2, length(colnames(dt_data)), 1),
                m = 2
            )
            
            ## for each combination, do a stat test
            for (combination_i in seq(1, ncol(combinations), 1)) {
                # print(combination_i)

                ## retreive the cluster names
                cluster_name1 <- colnames(dt_data)[combinations[1,combination_i]]
                cluster_name2 <- colnames(dt_data)[combinations[2,combination_i]]

                ## create the contingency table, removing the unknown data
                contingency_table <- dt_data[
                    variable != "unknown", 
                    c(
                        cluster_name1,
                        cluster_name2
                    ),
                    with = F
                ]


                ## transform to matrix
                contingency_table <- as.matrix(
                    contingency_table, 

                )

                # print(i)
                # print(contingency_table)

                ## add the row names
                rownames(contingency_table) <- unlist(dt_data[variable != "unknown", "variable"])



                if (nrow(contingency_table) > 2) {
                    contingency_table <- contingency_table[rownames(contingency_table) != 'deletion',]
                }                



                print(contingency_table)

                ## performe the statistical test
                statitiscal_results <- fisher.test(
                    x = contingency_table,
                    simulate.p.value = TRUE
                )
                #statitiscal_results <- chisq.test(
                #    x = contingency_table
                #)

                print(statitiscal_results)
                # print(statitiscal_results[["p.value"]])


                ## extract the pvalue and the clusters names used for the stat test
                pvalue <- as.data.table(
                    list(
                        cluster_name1,
                        cluster_name2,
                        statitiscal_results[["p.value"]]
                    ),
                )[, setnames(
                    .SD,
                    c("V1", "V2", "V3"),
                    colnames(dt_pvalues)
                )]
                
                ## add the pvalue to the pvalues data table 
                dt_pvalues <- rbind(
                    dt_pvalues,
                    pvalue
                )
                
            }

            ## put the summary table into the output list
            output_list[[i]][["data"]] <- dt_data
            output_list[[i]][["p_values"]] <- dt_pvalues
        }

        
    }

    ## return the ouput
    return(output_list)
}


###################################################
###################################################


##### FORMATE AND MERGE ALL THE CLINICAL DATA

## from the cluster column, extract cluster name and entity information for each sample
#d_subclusters <- d_subclusters[
#    ,
#     `:=` ( 
#        'entity' = sapply(cluster, function(x) str_split(x, pattern = '_cluster_')[[1]][2])
#    ) 
#][
#    ,
#    !c('cluster')
#]


## extract clinical data from the variable to keep
d_clinical_filtered <- unique(
    merge(
        d_clinical[, !c("Project_ID", "File_Name"), with = F],
        d_clinical_external,
        by.x = "Case_ID",
        by.y = "Case",
        all.x = T
    )
)

## merge with purity data
d_clinical_filtered <- merge(
    d_clinical_filtered,
    d_purity,
    by.x = "Case_ID",
    by.y = "NAME",
    all.x = T
)

## merge with cluster data
d_clinical_filtered <- merge(
    d_clinical_filtered,
    d_clusters,
    by.x = "Case_ID",
    by.y = "sample_id",
    all.y = T
)

## merge with the entity data
#d_clinical_filtered <- merge(
#    d_clinical_filtered,
#    d_subclusters,
#    by.x = "Case_ID",
#    by.y = "sample_id",
#    all.x = T
#)

## remove the EGFR column from the d_alteration data
d_alteration[, EGFR := NULL]

### merge with the alteration data
d_clinical_filtered <- merge(
    d_clinical_filtered,
    d_alteration,
    by.x = "Case_ID",
    by.y = "sample_id",
    all.x = T
)

## create a variable called "classification" describing the glioma classification : IDH-mut codel, IDH-mut noncodel and IDHwt
d_clinical_filtered[`IDH status` == "WT", "classification" := "IDHwt"]
d_clinical_filtered[(`IDH status` == "Mutant") & (`1p/19q codeletion` == "codel"), "classification" := "IDHmut_CODEL"]
d_clinical_filtered[(`IDH status` == "Mutant") & (`1p/19q codeletion` == "non-codel"), "classification" := "IDHmut_nonCODEL"]


###################################################
###################################################


## extract variables which are qualitatives
qualitative_variables <- c(
    "gender",
    "Histology",
    "Grade",
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
    "MYBL1"
)

#gaba_gluta_genes <- c(
#    'GNB1','GABRD','ALDH4A1','GRIK3','RIMKLA','SLC1A7','GNG12','PRKACB','GNG5','GNAI3','GLUL','PLA2G4A','GNG4','ADCY3','CAD','PPP3R1','GFPT1','KCNJ3','GAD1','GLS','PLCL1','TRAK2','CPS1','ITPR1','GRM7','SLC6A11','SLC6A1','SLC38A3','GNAI2','GRM2','CACNA1D','NIT2','ADCY5','TRPC1','PLD1','GNB4','NAT8L','GABRG1','GABRA2','GABRA4','GABRB1','PPAT','PPP3CA','GRIA2','ADCY2','SLC1A3','HOMER1','GRIA1','GABRB2','GABRA1','GABRG2','GFPT2','ALDH5A1','GABBR1','ITPR3','GRM4','GRIK2','DDO','GRM1','ADCY1','ASL','GNAI1','GRM3','GNG11','ASNS','GNB2','GRM8','PPP3CC','ADCY8','GPT','SLC1A1','GNAQ','GABBR2','GRIN3A','GNG10','ASS1','GRIN1','CACNA1B','GAD2','PPP3CB','GLUD1','GOT1','SLC17A6','SLC1A2','FOLH1','ASRGL1','GNG3','PLCB3','GRK2','SHANK2','GRM5','GRIA4','GRIK4','SLC6A13','SLC6A12','CACNA1C','GNB3','RIMKLB','GABARAPL1','GRIN2B','ITPR2','SLC38A1','SLC38A2','ADCY6','SLC17A8','ADCY4','GNG2','GPHN','GABRB3','GABRA5','GABRG3','PLCB2','JMJD7-PLA2G4B','GNB5','HOMER2','ADCY9','ABAT','GRIN2A','PRKCB','MAPK3','GPT2','ADCY7','GNAO1','GOT2','GABARAPL2','ASPA','PLD2','DLG4','GABARAP','HAP1','NSF','PRKCA','GRIN2C','DLGAP1','GNG7','CACNA1A','PRKACA','SLC1A6','HOMER3','GRIK5','PLA2G4C','GRIN2D','SLC17A7','IL4I1','SHANK1','PRKCG','PLCB1','PLCB4','SRC','SLC32A1','SLC12A5','GNAS','GRIK1','KCNJ6','MAPK1','GRK3','ADSL','SLC38A5','GLUD2','GRIA3','GABRE','GABRA3','GABRQ','MYB','EGFR','FGFR','MYBL1','CDKN2A','CDKN2B'
#)

#qualitative_variables <- c(
#    qualitative_variables,
#    gaba_gluta_genes
#)


## what variable are quantitative variable
quantitative_variables <- c(
    "age_at_index",
    "TumorPurity"
)

## extract data for each clusters
data_list <- extract_data(
    CLUSTER_GROUP_INPUT = d_clusters,
    CLINICAL_DATA_INPUT = d_clinical_filtered,
    QUALITATIVE_VARIABLE_INPUT = qualitative_variables,
    QUANTITATIVE_VARIABLE_INPUT = quantitative_variables
)
## statistical test
final_data <- statistic(
    EPIDEMIO_DATA_LIST_INPUT = data_list,
    QUALITATIVE_VARIABLE_INPUT = qualitative_variables,
    QUANTITATIVE_VARIABLE_INPUT = quantitative_variables)


###################################################
###################################################


##### WRITE THE OUTPUT


## extract tbe data of interest associated with all the samples
d_clinical_filtered <- d_clinical_filtered[, c("Case_ID","cluster", qualitative_variables, quantitative_variables), with = F]

## write the data table that contain all the epidemiological data
fwrite(
    x = d_clinical_filtered,
    file = snakemake@output[["epidemiological_data_all"]],
    sep = ","
)

## function for writting the output
write_data_list <- function(
    DATA_LIST_INPUT, # data list generated with the statistic function
    DIRECTORY_PATH_INPUT
) {

    ## extract the variable names
    variables <- names(DATA_LIST_INPUT)

    ## for each data of the list
    for (i_variable in seq(1, length(variables), 1)) {

        ## file name
        variable_name_for_file <- str_replace_all(
            variables[i_variable],
            pattern = "[ ]|[/]|[-]",
            replacement = "_"
        )

        # variable_name_for_file <- str_replace_all(
        #     variables[i_variable],
        #     pattern = "[/]",
        #     replacement = "_"
        # )

        ## file path for saving
        path_file_data <- paste(
            DIRECTORY_PATH_INPUT,
            "/",
            variable_name_for_file,
            "_table.csv",
            sep = ""
        )

        path_file_pvalues <- paste(
            DIRECTORY_PATH_INPUT,
            "/",
            variable_name_for_file,
            "_pvalues.csv",
            sep = ""
        )

        # write the data table
        fwrite(
            x = DATA_LIST_INPUT[[variables[i_variable]]][["data"]],
            file = path_file_data,
            sep = ","
        )

        # write the pvalues table
        fwrite(
            x = DATA_LIST_INPUT[[variables[i_variable]]][["p_values"]],
            file = path_file_pvalues,
            sep = ","
        )
    }
}

## create the directory path that will contain all the data and p values
data_path <- paste(
    snakemake@output[["epidemio_dir"]],
    "data",
    sep = "/"
)

## create the output directory from the data path
dir.create(
    data_path,
    recursive = T
)

## write the p value data into the output directory
write_data_list(
    DATA_LIST_INPUT = final_data,
    DIRECTORY_PATH_INPUT = data_path
)



###################################################
###################################################


##### MERGE DATA INTO A SAME DATA TABLE

merge_datalist <- function(
    DATA_LIST_INPUT,
    CLUSTER_COUNTING_INPUT,
    DIRECTORY_PATH_INPUT
) {

    ## extract the variable names
    variables <- names(DATA_LIST_INPUT)

    ## data output
    output <- c()

    ## for each data of the list
    for (i_variable in seq(1, length(variables), 1)) {

        ## file name
        variable_name_for_file <- str_replace_all(
            variables[i_variable],
            pattern = "[ ]|[/]|[-]",
            replacement = "_"
        )

        ## extract the data associated with the variable
        variable_data <- DATA_LIST_INPUT[[variables[i_variable]]][["data"]]

        ## put the header on the top of the ouput file and add the number total
        if (i_variable == 1) {
            header_output <- colnames(variable_data)

            ## create the header with the percentage of sample for each cluster
            header_output <- sapply(
                X = header_output,
                function(x) {
                    if (x == "variable"){
                        return("category")
                    }
                    else {
                        output_name <- paste(
                            x,
                            " (n = ",
                            toString(CLUSTER_COUNTING_INPUT[x]),
                            ")",
                            sep = ""
                        )
                        return(output_name)
                    }

                }
            )

            ## add "variable" in the begining of the header
            header_output <- c(
                "variable",
                header_output
            )

            ## add "," between each variable name
            header_output <- paste(
                header_output,
                collapse = ","
            )

            ### add a "," at the begining
            #header_output <- paste(
            #    ",",
            #    header_output,
            #    sep = ""
            #)

            ## add to the output
            output <- c(output, header_output)

        }

        ## put the variable name into the output 
        output <- c(output, variable_name_for_file)

        ## output each line of the data
        for (i_data_line in seq(1, nrow(variable_data), 1)) {
            data_line <- unlist(variable_data[i_data_line,])

            ## create the new line from the data_line
            new_data_line <- ""
            for (i_cluster_name in seq(1, length(names(data_line)), 1)) {

                ## retreive the cluster names
                cluster_name <- names(data_line)[i_cluster_name]

                ## if it is the variable names
                if (cluster_name %in% c("variable", "summary")) {
                    cluster_name_value <- as.character(data_line[names(data_line) == cluster_name])

                    ## save the nature of the variable, if it is a variable or a summary
                    nature <- cluster_name
                }

                else {

                    if (nature == "summary") {
                        ## retreive the cluster value
                        cluster_name_value <- as.numeric(data_line[names(data_line) == cluster_name])

                        cluster_name_value <- format(round(cluster_name_value, 2), nsmall = 2)
                    }

                    ## if == variable, calculate the percentage for each clusters
                    if (nature == "variable") {
                        ## retreive the cluster value
                        cluster_name_value <- as.numeric(data_line[names(data_line) == cluster_name])

                        ## calculate the percentage
                        percentage <- cluster_name_value/CLUSTER_COUNTING_INPUT[[cluster_name]] * 100
                        percentage <- format(round(percentage, 1), nsmall = 1)

                        cluster_name_value <- paste(
                            cluster_name_value,
                            " (",
                            percentage,
                            "%)",
                            sep = ""
                        )

                    }
                   
                }
                
                ## collapse all the value to the new data line
                new_data_line <- paste(
                    new_data_line,
                    cluster_name_value,
                    sep = ","
                )
            }
            ## add into the output
            output <- c(output, new_data_line)
        }
    }

    ### add some empty string in the beggining for the column
    #output <- c(
    #    paste(paste("V", seq(1, 6, 1), sep = ""), collapse = ","),
    #    output
    #)

    ## add a "\n" between each element of the output
    output <- paste(
        output,
        collapse = "\n"
    )


    ## write the output
    cat(
        output,
        file = DIRECTORY_PATH_INPUT
    )
}




## merge all the data and add the percentage
merge_datalist(
    DATA_LIST_INPUT = final_data,
    CLUSTER_COUNTING_INPUT = cluster_counting,
    DIRECTORY_PATH_INPUT = snakemake@output[["epidemiological_data_formated"]]
)



###############################
## for gender and age data
###############################

## extract the data that correspond to the gender, sex and age
epidemio_data_of_interest <- final_data[names(final_data) %in% c("gender", "age_at_index")]

## merge all the data and add the percentage
merge_datalist(
    DATA_LIST_INPUT = epidemio_data_of_interest,
    CLUSTER_COUNTING_INPUT = cluster_counting,
    DIRECTORY_PATH_INPUT = snakemake@output[["epidemiological_data_formated_subset"]]
)


