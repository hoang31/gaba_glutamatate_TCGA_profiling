

#########################################################
#########################################################


##### CALCULATE THE SURVIVAL FOR EACH GLIOMA CLUSTERS


#########################################################
#########################################################


##### LOAD THE LIBRARIES
library(data.table)
library(survival)
library(survminer)
library(RColorBrewer)
library(stringr)

#########################################################
#########################################################


##### LOAD DATA

### load the cluster data
cluster_group <- fread(
    file =  snakemake@input[["cluster_group"]],
    sep = ",",
)

## load the subgroup data
epidemio_dt <- fread(
    snakemake@input[['subgroup_epidemiological_data_all']],
    sep = ','
)

## vector of the cluster names
cluster_vector <- unique(unlist(epidemio_dt[, cluster]))

## vector of the subcluster names
subcluster_vector <- unique(unlist(epidemio_dt[, classification2]))

## load the clinical data
clinical_data <- fread(
    file = snakemake@input[["clinical_data"]],
    sep = ","
)[, !c("File_Name", "Project_ID"), with = F]

## load the external clinical data
d_clinical_TCGA_external <- fread(
    input = snakemake@input[["TCGA_clinical_EXTERNAL"]],
    sep = ",",
    header = T
)

## recuperation of the follow up data
follow_up_data <- fread(
    file = snakemake@input[["clinical_data_full"]],
    sep = "\t"
)

variables_to_keep <- c(
    "case_submitter_id", # corresponding to the primary key
    "days_to_last_follow_up"
)

## keep follow up data
follow_up_data <- unique(follow_up_data[, variables_to_keep, with = F])

## remove duplicate rows
clinical_data <- unique(clinical_data)

## extract the wild card from the snakemake rule
wcard <- snakemake@params[["wcard"]]

## load survival function from the utils.R scripts
source(snakemake@params[["utils_R"]])

## create the directory
dir.create(snakemake@output[["survival_curve_dir"]])

#print('######################')
#print(table(cluster_group[, cluster]))

#########################################################
#########################################################

##### Cluster data


##### FORMATE THE DATA for the survival analysis

## variable to keep
relevant_informations <- c(
    "sample_id", 
    "cluster",
    "vital_status",
    "days_to_death"
)

## merge the cluster information and the clinical data and extract the relevant variables for the survival analysis
survival_data <- merge(
    x = cluster_group,
    y = clinical_data,
    by.x = "sample_id",
    by.y = "Case_ID", 
    all.x = T
)[, relevant_informations, with = F]

## merge with the follow up data
survival_data <- merge(
    x = survival_data,
    y = follow_up_data,
    by.x = "sample_id",
    by.y = "case_submitter_id", 
    all.x = T
)

## replace the na value of the days_to_death column by the days to last follow up value for the alive persons and then remove the days_to_last_follow_up column
survival_data[vital_status == "Alive", days_to_death := days_to_last_follow_up ][,days_to_last_follow_up := NULL]

## remove the rows with no vital status
survival_data <- survival_data[!(vital_status == ""), ]

## remove the rows with no vital status
survival_data <- survival_data[!(vital_status == ""), ]

## remove the rows with no days_to_death
survival_data <- survival_data[!(days_to_death == "--"), ]

## transformorm the "vital_status" column as binary variable
survival_data[, vital_status := ifelse(vital_status == "Dead", 1, 0)]

## days_to_death : transform days to months
survival_data[, days_to_death := sapply(
    days_to_death, 
    function(x) {
        return(as.numeric(x) * 12 / 365)
    }
)]

## rename the column names
survival_data <- survival_data[, setnames(
    .SD, 
    c("days_to_death", "vital_status"), 
    c("time", "status"))]


##### SURVIVAL ANLYSIS

## colors that will be used for the generation of the plots
color_palette <- brewer.pal(9, "Set1")

## create a datatable containing the association of colors for each cluster name
color_data <- data.table(
    cluster = sort(cluster_vector),
    color = color_palette[1:(length(cluster_vector))]
)

## create the survival curve calling the function
create_survival_curve(
    FIT_MODEL_INPUT = survfit(Surv(time, status) ~ cluster, data = survival_data),
    FILE_PATH_INPUT = paste(
        snakemake@output[["survival_curve_dir"]],
        "survival_curve.svg",
        sep = "/"
    ),
    COLOR_PALETTE_INPUT = unlist(color_data[,color]),
    DISPLAY_PVALUE = F,
    WIDTH_INPUT = 12,
    HEIGTH_INPUT = 12
)


##########
## Generate the p values for each cluster comparison
##########

## retreive all the cluster name into a vector
cluster_name_vector <- unique(unlist(cluster_group[, cluster]))

## initialize the data table that will contain the pvalues for each combination
pval_all_dt <- data.table(
    "cluster1" = character(),
    "cluster2" = character(),
    "pval" = numeric()
)

## do the combination of two cluster
cluster_combination <- combn(
    cluster_name_vector,
    m = 2
)

## for each combination, extract the data and generate the survival data
for (i_combination in seq(1, ncol(cluster_combination))) {

    ## retreive the cluster 1
    cluster_name1 <- cluster_combination[1, i_combination]
    cluster_name2 <- cluster_combination[2, i_combination]

    ## copy the subset data of the survival data
    subdata_survival <- copy(survival_data)

    ## extract data related to the combination cluster
    subdata_survival <- subdata_survival[cluster %in% c(cluster_name1, cluster_name2),]

    ## generate the survival data
    sub_fitdata <- survfit(Surv(time, status) ~ cluster, data = subdata_survival)

    ## retreive the pvalues
    pval <- surv_pvalue(sub_fitdata)$pval

    ## put the cluster combination names associated with their pvalue in a data table
    pval_dt <- data.table(
        "cluster1" = cluster_name1,
        "cluster2" = cluster_name2,
        "pval" = pval
    )

    ## rbind the pval_dt with the data table that contain all the pvalue
    pval_all_dt <- rbind(
        pval_all_dt,
        pval_dt
    )
}

## write the data
fwrite(
    pval_all_dt,
    paste(
        snakemake@output[["survival_curve_dir"]],
        "/",
        "pvalues.csv",
        sep = ""
    ),
    sep = ","
)

#########################################################
#########################################################

##### SubCluster data

#print(survival_data)
#print(epidemio_dt)

#epidemio_dt <- epidemio_dt[subcluster != 'IDHmut_nonCODEL_cluster_IDHwt_astrocytoma',]

#print(sort(unique(unlist(epidemio_dt[, subcluster]))))
##break




for (i_subcluster in seq(1, length(subcluster_vector), 1)) {

    ## extract the subcluster name
    subcluster_name <- subcluster_vector[i_subcluster]

    ## extract the sample id associated with the subcluster name
    sample_id_vector <- epidemio_dt[classification2 == subcluster_name, Case_ID]

    ## extract the survival data associated with the sample id vector
    survival_subdata <- copy(
        survival_data[sample_id %in% sample_id_vector, ]
    )

    ## generate the surival curve
    create_survival_curve(
        FIT_MODEL_INPUT = survfit(Surv(time, status) ~ cluster, data = survival_subdata),
        FILE_PATH_INPUT = paste(
            snakemake@output[["survival_curve_dir"]],
            "/survival_curve_",
            subcluster_name,
            ".png",
            sep = ""
        ),
        COLOR_PALETTE_INPUT = color_palette[1:5],
        DISPLAY_PVALUE = F,
        WIDTH_INPUT = 1400,
        HEIGTH_INPUT = 800
    )
    #print(subcluster_name)
    #print(table(survival_subdata[, cluster]))
    #print(survival_subdata)
}


#########################################################
#########################################################

## Do a cox regression analysis corrected on the age, the sex and the histology

## merge the survival data with the variables that will be taking in account for the regression correction

## load the epidemiological_data
epidemiological_data <- fread(
    snakemake@input[['epidemiological_data_all']],
    sep = ','
)

## put the variable to keep
variable_to_keep <- c(
    'Case_ID',
    'Histology',
    'age_at_index',
    'Grade',
    "gender"
    #"classification",
)

## merge the survival data with the variable to keep
survival_data_cox <- merge(
    survival_data,
    epidemiological_data[, variable_to_keep, with = F],
    by.x = 'sample_id',
    by.y = 'Case_ID',
    all.x = T
)

## rename the blank values to unknown values
survival_data_cox[survival_data_cox == ""] <- "unknown"



## merge the survival data with the variable to keep
#survival_data_cox <- merge(
#    survival_data,
#    epidemiological_data,
#    by.x = 'sample_id',
#    by.y = 'Case_ID',
#    all.x = T
#)

#fwrite(
#    x = survival_data_cox,
#    file = "data/survival_data_epidemio.csv",
#    sep = ","
#)

#Sys.sleep(100000)
#break


## extract only the IDHwt and MIXED gliomas
#survival_data_cox <- survival_data_cox[cluster %in% c('IDHwt', 'MIXED')]
#survival_data_cox <- survival_data_cox[cluster %in% c('NT-1', 'NT-2')]
#survival_data_cox <- survival_data_cox[cluster %in% c('NT-4', 'NT-2')]


## function for getting the results from the cox regression model
get_cox_regression_results <- function(
    cox_model
) {

    ## do the sumarize
    cox_model_summary  <- summary(cox_model)

    ## extract the coefficient
    values <- as.data.table(cox_model_summary$coefficients, keep.rownames = 'name')
    name <- unlist(values[, 'name', with = F])
    beta <- format(round(unlist(values[, 'coef', with = F]), 3), nsmall = 3)
    hazard_ratio <- format(round(unlist(values[, 'exp(coef)', with = F]), 3), nsmall = 3)
    pvalue <- unlist(values[, 'Pr(>|z|)', with = F])


    ### extract confidence intervals
    interval_confidence <- as.data.table(cox_model_summary$conf.int, keep.rownames = 'name')
    confidence_interval_lower <- format(round(unlist(interval_confidence[, 'lower .95', with = F]), 3), nsmall = 3)
    confidence_interval_upper <- format(round(unlist(interval_confidence[, 'upper .95', with = F]), 3), nsmall = 3) 

    ## formate the interval confidence
    interval_confidence <- paste(
        hazard_ratio,
        ' (',
        confidence_interval_lower,
        ' - ',
        confidence_interval_upper,
        ')',
        sep = ''
    )
    
    ## put together all the data into a data table
    res <- data.table(
        'name' = name,
        'beta' = beta,
        'HR (95% CI for HR)' = interval_confidence,
        'pvalue' = pvalue
    )

    ## return the output
    return(res)
} 

## covariate of interest
covariates <- c(
    'cluster',
    'gender',
    'age_at_index',
    'Grade',
    'cluster + gender + age_at_index + Grade'
    #"classification",
    #'cluster + gender + age_at_index + Grade + classification'

)

## generate formulas for each covariate of interest
cox_formulas <- sapply(
    covariates,
    function(x) 
        as.formula(paste('Surv(time, status) ~', x))
)

                        
## generate the cox regression model for each formula
cox_models <- lapply(
    cox_formulas,
    function(x)
        {
            coxph(x, data = survival_data_cox)
        }
    )



## extract the data and the metrics for each regression model
results <- lapply(
    cox_models,
    get_cox_regression_results
)

####################################################
## Write the results into files
####################################################

for (i_results in seq(1, length(results), 1)) {
    
    ## extract names
    names <- names(results[i_results])

    ## generate the file path for the file
    file_path <- paste(
        snakemake@output[['survival_curve_dir']],
        '/cox_regression ',
        names,
        '.csv',
        sep = ''
    )

    output <- data.table(results[i_results][[1]])
 
    ## write the file into a file 
    fwrite(
        x = output,
        file = file_path,
        sep = ','
    )
}


