
#########################################################
#########################################################


##### SCRIPT FOR CLINICAL DIFFERENCES BETWEEN EACH DISTINCT ENTITY FOR EACH CLUSTER


#########################################################
#########################################################


##### LOAD PACKAGES
library(data.table)
library(tidyverse)


#########################################################
#########################################################

##### LOAD THE DATA

## load the epidemiological data
epidemio_data <- fread(
    file = snakemake@input[['epidemiological_data_all']]
)

## load external scripts
source(snakemake@params[["utils_R"]])

#########################################################
## Extract each entity for each cluster
#########################################################

## initialize the classification2 variable 
epidemio_data[
    ,
    classification2 := 'unknown'
]

## IDHwt astrocytomas 
epidemio_data[
    ((Histology == 'glioblastoma') & (`IDH status` == 'WT')),
    classification2 := 'IDHwt_glioblastoma'
]
epidemio_data[
    ((`IDH status` == 'WT') & ((`TERT promoter status` == 'Mutant') | (EGFR == 'amplification') | (`Chr 7 gain/Chr 10 loss` == 'Gain chr 7 & loss chr 10'))),
    classification2 := 'IDHwt_glioblastoma'
]

epidemio_data[
    ((`classification2` == 'unknown') & (`IDH status` == 'WT')),
    classification2 := 'IDHwt_astrocytoma'
]

epidemio_data[
    ((`IDH status` == 'WT') & (classification2 == 'unknown')),
    classification2 := 'IDHwt_unknown_astrocytoma'
]

## IDH mutated astrocytomas
epidemio_data[
    ((`IDH status` == 'Mutant') & (Grade == 'G2')),
    classification2 := 'IDHmut_difuse_astrocytoma'
]
epidemio_data[
    ((`IDH status` == 'Mutant') & (Grade == 'G3')),
    classification2 := 'IDHmut_anaplasic_astrocytoma'
]
epidemio_data[
    ((`IDH status` == 'Mutant') & (Grade == 'G4')),
    classification2 := 'IDHmut_glioblastoma'
]
epidemio_data[
    ((`IDH status` == 'Mutant') & (classification2 == 'unknown')),
    classification2 := 'IDHmut_unknown_astrocytoma'
]

## oligodenndroglioma
epidemio_data[
    ((`1p/19q codeletion` == 'codel') & (`IDH status` == 'Mutant')),
    classification2 := 'oligodendroglioma'
]





#variable <- c(
#    'IDH status',
#    'Grade',
#    'Histology',
#    'classification2'
#)
#test <- epidemio_data[, variable, with = F][order(`IDH status`, `Grade`),]
#setDF(test)
#print(test)
#break






## merge the classification2 variable with the cluster name
epidemio_data[
    ,
    subcluster := paste(
        cluster,
        '_cluster_',
        classification2,
        sep = ''
    )
]

## extract the subcluster data and rename the column
d_subclusters <- epidemio_data[,c('Case_ID', 'subcluster')][, setnames(.SD, c('Case_ID', 'subcluster'), c('sample_id', 'cluster'))]

cat("\n\n##################################################")
cluster_counting <- table(d_subclusters[,"cluster"])
print(cluster_counting)
subcluster_counting <- table(d_subclusters[,"cluster"])
print(subcluster_counting)
cat("###################################################\n\n")


#########################################################
#########################################################


##### FUNCTIONS

## function for extracting the id for each entity for each cluster group
extract_sample_id <- function(
    cluster_data_input # contains all the cluster informations and the sample ids
) {

    ## copy the cluster data table
    cluster_data_input <- copy(cluster_data_input)

    ## identify the cluster and histology names for each sample
    cluster_data_input[, cluster_name := sapply(cluster, function(x) str_split(x, pattern = "_cluster_")[[1]][1])]
    cluster_data_input[, histology_name := sapply(cluster, function(x) str_split(x, pattern = "_cluster_")[[1]][2])]

    ## retreive the unique vector of cluster and histology names
    cluster_name_vector <- unique(unlist(cluster_data_input[, "cluster_name"]))
    histology_name_vector <- unique(unlist(cluster_data_input[, "histology_name"]))

    ## initialize the list that will contain the sample id of each cluster and histology
    sample_id_list <- list()

    ## for each cluster name and histology name
    for (i_cluster_name in seq(1, length(cluster_name_vector), 1)) {
        
        ## extract the cluster name at the i-th position
        cluster_name2 <- cluster_name_vector[i_cluster_name]

        for (i_histology_name in seq(1, length(histology_name_vector), 1)) {

            ## extract the histology name at the i-th position
            histology_name2 <- histology_name_vector[i_histology_name]
            
            # cat("----------------------",cluster_name2, histology_name2, "\n")
            ## extract the sample id correcponding to the cluster name and the histology names
            sample_id_vector <- unlist(cluster_data_input[(cluster_name == cluster_name2) & (histology_name == histology_name2), "sample_id" ])

            ## remove the names
            names(sample_id_vector) <- NULL

            # print(sample_id_vector)
            ## add sample id in the list
            if (length(sample_id_vector) > 0) {
                sample_id_list[[histology_name2]][[cluster_name2]] <- sample_id_vector
            }
            else {
                sample_id_list[[histology_name2]][[cluster_name2]] <- NA
            }
        }
    }

    ## return the sample_id list
    return(sample_id_list)

}


## function for extracing informations associated with the sample_id list
extract_informations <- function(
    sample_id_list_input, # list of sample id generated by the "extract_sample_id" function
    cluster_data_input, # contains all the cluster informations and the sample ids
    clinical_data_input, # clinical data table
    qualitative_variable_input,
    quantitative_variable_input
) {

    ## copy the cluster data table
    cluster_data_input <- copy(cluster_data_input)

    ## identify the cluster and histology names for each sample
    cluster_data_input[, cluster_name := sapply(cluster, function(x) str_split(x, pattern = "_cluster_")[[1]][1])]
    cluster_data_input[, histology_name := sapply(cluster, function(x) str_split(x, pattern = "_cluster_")[[1]][2])]

    ## retreive the unique vector of cluster and histology names
    cluster_name_vector <- unique(unlist(cluster_data_input[, "cluster_name"]))
    histology_name_vector <- unique(unlist(cluster_data_input[, "histology_name"]))

    ## merge the qualitative and quantitative variables into a same vector
    variable_vector <- c(
        qualitative_variable_input,
        quantitative_variable_input
    )

    ## initialize the listthat will contain the output
    data_list <- list()

    ## for each variable
    for (i_variable_vector in seq(1, length(variable_vector), 1)) {
        
        ## extract the variable name
        variable_name <- variable_vector[i_variable_vector]

        ## extract all the data associated with the variable
        variable_data <- clinical_data_input[,c("Case_ID", variable_name), with = F]

        if (variable_name %in% qualitative_variable_input) {
            
            ## remove NA and empty category and replace them by 'unknow'
            variable_data[is.na(variable_data)] <- "unknown"
            variable_data[variable_data == ""] <- "unknown"

            ## extract all the categories associated with the variable
            category_vector <- unique(unlist(variable_data[, variable_name, with = F]))

            ## initialize a data table containing all the category vector
            category_dt <- data.table(category = category_vector)
        }

        if (variable_name %in% quantitative_variable_input) {

            ## do the summary and initialize the data table that will contain the variable category (mean, min, median, max, etc)) 
            numeric_vector <- summary(as.numeric(unlist(variable_data[, variable_name, with = F])))
            numeric_vector <- as.data.table(as.matrix(numeric_vector), keep.rownames = T)[, setnames(.SD, c("rn", "V1"), c("category", "N"))]
            category_dt <- numeric_vector[, !c("N"), with = F]
        }

        #cat("\n")
        #print("#############################################################################")
        #print(variable_name)
        #print(category_dt)
        #next

        ## for each cluster name and histology name
        for (i_cluster_name in seq(1, length(cluster_name_vector), 1)) {
            
            ## extract the cluster name at the i-th position
            cluster_name2 <- cluster_name_vector[i_cluster_name]
            
            for (i_histology_name in seq(1, length(histology_name_vector), 1)) {

                ## extract the histology name at the i-th position
                histology_name2 <- histology_name_vector[i_histology_name]
                            
                #cat("------------------------",cluster_name2, histology_name2, variable_name, "\n")

                ## extract the sample id for each cluster and histology names
                sub_cluster_sample_id <- sample_id_list_input[[histology_name2]][[cluster_name2]]
                
                ## for each variable, extract information associated with the variable
                sub_cluster_variable_data <- variable_data[Case_ID %in% sub_cluster_sample_id, variable_name, with = F]

                #print(sub_cluster_variable_data)

                ## if the variable is qualitative or quantitative
                if (variable_name %in% qualitative_variable_input) {
                    ## replace the na by unknown
                    sub_cluster_variable_data[is.na(sub_cluster_variable_data)] <- "unknown"

                    ## if there are not samples
                    if (nrow(sub_cluster_variable_data) == 0) {
                        ## create a empty data table
                        sub_cluster_variable_data_table <- as.data.table(NULL)
                    }

                    ## if there are samples
                    else {
                        ## do the count of each categories for the variable name
                        sub_cluster_variable_data_table <- as.data.table(table(sub_cluster_variable_data))
                    }
                }

                if (variable_name %in% quantitative_variable_input) {

                    ## tranform to numeric
                    sub_cluster_variable_data <- as.numeric(unlist(sub_cluster_variable_data))

                    ## do the summary and transform the output to a data table
                    sub_cluster_variable_data_table <- as.data.table(as.matrix(summary(sub_cluster_variable_data)), keep.rownames = T)[, setnames(.SD, c("rn", "V1"), c("sub_cluster_variable_data", "N"))]
                }

                ## if there are sample of the group
                if (nrow(sub_cluster_variable_data_table) > 0) {
                    ## merge the count table with the category datatable
                    sub_cluster_variable_data_table <- merge(
                        category_dt,
                        sub_cluster_variable_data_table,
                        by.x = "category",
                        by.y = "sub_cluster_variable_data",
                        all.x = T
                    )

                    ## put 0 to the NA values
                    sub_cluster_variable_data_table[is.na(sub_cluster_variable_data_table)] <- 0
                }
                ## if there are not samples 
                else {
                    ## copy and put 0 to all the values
                    sub_cluster_variable_data_table <- copy(category_dt)
                    sub_cluster_variable_data_table[, N := 0]
                }

                ## rename the column by the cluster type
                sub_cluster_variable_data_table <- sub_cluster_variable_data_table[, setnames(.SD, "N", cluster_name2)]
                sub_cluster_variable_data_table <- sub_cluster_variable_data_table[, setnames(.SD, "category", variable_name)]

                ## add the count table in the data_list output
                data_list[[histology_name2]][[variable_name]][[cluster_name2]] <- sub_cluster_variable_data_table

            }
        }
    }

    ## return the output
    return(data_list)
}

#########################################################
#########################################################


##### Extrac the information of interest

## create a variable called "classification" describing the glioma classification : IDH-mut codel, IDH-mut noncodel and IDHwt
epidemio_data[`IDH status` == "WT", "classification" := "IDHwt"]
epidemio_data[(`IDH status` == "Mutant") & (`1p/19q codeletion` == "codel"), "classification" := "IDHmut_CODEL"]
epidemio_data[(`IDH status` == "Mutant") & (`1p/19q codeletion` == "non-codel"), "classification" := "IDHmut_nonCODEL"]

## replace the zero by wt
epidemio_data[epidemio_data == 0] <- 'wt'

## store the qualitative variables into a vector
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
    "MYBL1"
)


#gaba_gluta_genes <- c(
#    'GNB1','GABRD','ALDH4A1','GRIK3','RIMKLA','SLC1A7','GNG12','PRKACB','GNG5','GNAI3','GLUL','PLA2G4A','GNG4','ADCY3','CAD','PPP3R1','GFPT1','KCNJ3','GAD1','GLS','PLCL1','TRAK2','CPS1','ITPR1','GRM7','SLC6A11','SLC6A1','SLC38A3','GNAI2','GRM2','CACNA1D','NIT2','ADCY5','TRPC1','PLD1','GNB4','NAT8L','GABRG1','GABRA2','GABRA4','GABRB1','PPAT','PPP3CA','GRIA2','ADCY2','SLC1A3','HOMER1','GRIA1','GABRB2','GABRA1','GABRG2','GFPT2','ALDH5A1','GABBR1','ITPR3','GRM4','GRIK2','DDO','GRM1','ADCY1','ASL','GNAI1','GRM3','GNG11','ASNS','GNB2','GRM8','PPP3CC','ADCY8','GPT','SLC1A1','GNAQ','GABBR2','GRIN3A','GNG10','ASS1','GRIN1','CACNA1B','GAD2','PPP3CB','GLUD1','GOT1','SLC17A6','SLC1A2','FOLH1','ASRGL1','GNG3','PLCB3','GRK2','SHANK2','GRM5','GRIA4','GRIK4','SLC6A13','SLC6A12','CACNA1C','GNB3','RIMKLB','GABARAPL1','GRIN2B','ITPR2','SLC38A1','SLC38A2','ADCY6','SLC17A8','ADCY4','GNG2','GPHN','GABRB3','GABRA5','GABRG3','PLCB2','JMJD7-PLA2G4B','GNB5','HOMER2','ADCY9','ABAT','GRIN2A','PRKCB','MAPK3','GPT2','ADCY7','GNAO1','GOT2','GABARAPL2','ASPA','PLD2','DLG4','GABARAP','HAP1','NSF','PRKCA','GRIN2C','DLGAP1','GNG7','CACNA1A','PRKACA','SLC1A6','HOMER3','GRIK5','PLA2G4C','GRIN2D','SLC17A7','IL4I1','SHANK1','PRKCG','PLCB1','PLCB4','SRC','SLC32A1','SLC12A5','GNAS','GRIK1','KCNJ6','MAPK1','GRK3','ADSL','SLC38A5','GLUD2','GRIA3','GABRE','GABRA3','GABRQ','MYB','EGFR','FGFR','MYBL1','CDKN2A','CDKN2B'
#)

#qualitative_variables <- c(
#    qualitative_variables,
#    gaba_gluta_genes
#)


## store the quantitatives variables into a other vector
quantitative_variables <- c(
    "age_at_index",
    "TumorPurity"
)

## extract the sample ids of each subgroup
sample_id_list <- extract_sample_id(
    cluster_data_input = d_subclusters
)

## extract the information of each subgroup
clinical_data_list <- extract_informations(
    sample_id_list_input = sample_id_list, # list of sample id generated by the "extract_sample_id" function
    cluster_data_input = d_subclusters, # contains all the cluster informations and the sample ids
    clinical_data_input = epidemio_data, # clinical data table
    qualitative_variable_input = qualitative_variables,
    quantitative_variable_input = quantitative_variables
)

## extract the variable and histology names
histology_name_vector <-names(clinical_data_list)
variable_name_vector <- names(clinical_data_list[[1]])
cluster_name_vector <- names(cluster_counting)

## do the combination of all the histological names in pair
cluster_name_combination <- combn(
    cluster_name_vector,
    m = 2
)

output <- paste(
    ",,",
    paste(
        cluster_name_vector,
        collapse = ","
    ),
    sep = ""
)

## Initialize the output datatable list
output_dt_list <- list()

## for each variable
for (i_histology_name_vector in seq(1, length(histology_name_vector), 1)) {
    
    ## extract the histology name
    histology_name <- histology_name_vector[i_histology_name_vector]

    ## initialize the data table that will contain the output
    output_dt <- data.table(
        variable = NULL,
        MIXED = NULL,
        IDH_WT = NULL,
        IDH_MUT_noncodel = NULL,
        IDH_MUT_codel = NULL
    )

    cat("------------------------", histology_name, "------------------------", "\n")

    ## for the variable i and for each histology
    for (i_variable_name_vector in seq(1, length(variable_name_vector), 1))  {

        variable_name <- variable_name_vector[i_variable_name_vector]    
        
        ## merge all the sublist of each cluster associated with the same histoloy name
        data <- Reduce(
            x = clinical_data_list[[histology_name]][[variable_name]],
            function(d1, d2) {
                merge(
                    d1,
                    d2,
                    by = variable_name,
                    all.x = T
                )
            }
        )

        print(data)
        
        ################################
        ## to do : statistical test between each cluster
        ################################

        # ## do the statistical test for each combination of clusters
        #for (i_histology_combination in seq(1, ncol(cluster_name_combination))) {
        #    histology_combination <- cluster_name_combination[, i_histology_combination]

        #    print('-------------')
        #    print(variable_name)

        #    ## subset the data associated with the combination
        #    data_subset <- data[, histology_combination, with = F]

        # if (variable_name %in% qualitative_variables) {
        #     print('do the statistical test for qualitative variable')                
            
        #     print(data_subset)
        #     res <- chisq.test(data_subset)

        #     print(res)

        # }
            # if (variable_name %in% quantitative_variables) {
            #     print('do the statistical test for quantitative variable')                
            # }
        #}
    

        # next

        ## create a new row that contains the variable names
        new_row_variable <- colnames(data)
        new_row_variable[new_row_variable != variable_name] <- ""
        new_row_variable <- as.data.table(as.list(new_row_variable))

        ## add the new row containing the variable name into the data
        data2 <- as.data.table(
            rbindlist(
                list(
                    new_row_variable,
                    data
                ),
                fill = F
            )
        )

        # create a new row that will contain the variable name
        if (i_variable_name_vector == 1) {

            ## create a row containing the total
            row_total <- unlist(data[, lapply(.SD, sum), .SDcols = colnames(data)[colnames(data) != variable_name]])

            ## formate the total value
            row_total <- sapply(
                row_total,
                function(x) paste(
                    "(N = ",
                    x,
                    ")",
                    sep = ""
                )
            )
            row_total <- as.data.table(as.list(c(histology_name, row_total)))
            colnames(row_total)[colnames(row_total) == ""] <- "variable"

            data2 <- as.data.table(rbindlist(
                list(
                    row_total,
                    data2
                )
            ))
        }

        ## rename the names of the columns replacing the variable name by the string "variable"
        new_header <- colnames(data)[!(colnames(data) == variable_name)]
        new_header <- c("variable", new_header)
        colnames(data2) <- new_header
        
        ## merge the output data table with the new data
        output_dt <- rbind(
            output_dt,
            data2
        )
    }

    ## add the output_dt into a same list
    output_dt_list[[histology_name]] <- output_dt
}

## merge all the entity subgroup intp a same data table
output_merged <- Reduce(
    rbind, 
    output_dt_list
)

#################################################
## write the outputs
#################################################

## write the file into a file
fwrite(
    output_merged,
    snakemake@output[['subgroup_epidemiological_data_formated']],
    sep = ","
)

fwrite(
    epidemio_data,
    snakemake@output[['subgroup_epidemiological_data_all']],
    sep = ","
)
