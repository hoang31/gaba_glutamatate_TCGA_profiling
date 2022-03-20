
######################################
## Generate cleaned descriptive tables
######################################

######################################
## Load the libraries
######################################

library(data.table)
library(tidyverse)
library(sjPlot)


######################################
## create the output directory
######################################

## create the directory that will contain all the tables
dir.create(snakemake@output[["descriptive_table_dir"]])

######################################
## For the gender and age epidemiological data
######################################



### load the data into a data table
#data_dt <- fread(
#    "/home/hoang/Documents/projects/gaba_and_glutamate_pathways_in_glioma_v2/data/TCGA_IDHall_epidemiology/data/Karnofsky_Performance_Score_table.csv",,
#    sep = ","
    
#)[ variable != "unknown",]


#subset_data <- data_dt[, c("NT-4", "NT-2")]
#subset_data <- data_dt[, c("NT-3", "NT-2")]

#fisher.test(subset_data, simulate.p.value=TRUE)


#break



## load the epidemiological data
epidemio_data_dt <- fread(
    snakemake@input[["formated_epidemiological_data"]],
    sep = ",",
    header = T,
    fill = T
)

## sort the column
sorted_colnames <- c(
    "variable",
    #"category",
    sort(colnames(epidemio_data_dt)[!(colnames(epidemio_data_dt) %in% c("variable", "categorie"))])
)

## reanrange the order of the column
epidemio_data_dt <- epidemio_data_dt[,sorted_colnames, with = F]

## add a capital letter for the first letter
epidemio_data_dt <- epidemio_data_dt[, setnames(.SD, c("variable", "category"), c("Variable", "Category"))]

## save the table into a file in the directory
tab_df(
    epidemio_data_dt,
    file = paste(
        snakemake@output[["descriptive_table_dir"]],
        "/",
        "epidemio_table",
        ".html",
        sep = ""
    )
)

## save the table into a file in the directory
fwrite(
    epidemio_data_dt,
    paste(
        snakemake@output[["descriptive_table_dir"]],
        "/",
        "epidemio_table",
        ".csv",
        sep = ""
    ),
    sep = ","
)



######################################
## For the cox regression analysis
######################################

## load the path of all the files from the directory that contains the survival data
cox_file_path_vector <- list.files(
    snakemake@input[["survival_data_dir"]],
    full.names = T
)

## load only the name files
cox_file_name_vector <- list.files(
    snakemake@input[["survival_data_dir"]],
    full.names = F
)

## retreive only the path that contain cox regression data
cox_file_path_vector <- grep(
    x = cox_file_path_vector,
    pattern = "cox_regression",
    value = T
)

## retreive the name of the files related to cox data and reformate the name
cox_file_name_vector <- grep(
    x = cox_file_name_vector,
    pattern = "cox_regression",
    value = T
)
cox_file_name_vector <- str_replace(
    cox_file_name_vector,
    pattern = "cox_regression ",
    replacement = ""
)
cox_file_name_vector <- str_replace(
    cox_file_name_vector,
    pattern = ".csv",
    replacement = ""
)

## load all the cox data 
cox_data <- lapply(
    cox_file_path_vector,
    function(x) fread(
        x,
        sep = ","
    )
)

## rename the cox data
names(cox_data) <- cox_file_name_vector

## extract the data that contain the multivariate regression data, to do so, we will check if the data name is composed of " + " between covariate
multivariate_cox_file_name <- grep(
    cox_file_name_vector,
    pattern = "[+]",
    value = T
)

## extract all the covariates that compose the mutltivariate analysis
covariate_vector <- unlist(
    str_split(
        string = multivariate_cox_file_name,
        pattern = "[ + ]"
    )
)

## remove the blank values
covariate_vector <- covariate_vector[!(covariate_vector == "")]

## colname of the cox data
cox_data_colnames <- c(
    "covariate",
    "beta",
    "HR (95% CI for HR)",
    "pvalue"
)

## extract the multivariate data 
multivariate_data <- covariate_data <- as.data.table(cox_data[multivariate_cox_file_name])

## rename the multivariate data
colnames(multivariate_data) <- cox_data_colnames

## create the table that contain the univariate data in the same order of the covariate that compose the multivariate data
univariate_data <- data.table(
    "covariate" = character(),
    "beta"= character(),
    "HR (95% CI for HR)" = character(),
    "pvalue" = character()
)

for (i_covariate in seq(1, length(covariate_vector))) {

    ## retreive the covariate name
    covariate_name <- covariate_vector[i_covariate]

    ## retreive the data associated with the covariate name
    covariate_data <- as.data.table(cox_data[covariate_name])
    
    ## rename the data
    colnames(covariate_data) <- cox_data_colnames
    
    ## add the covariate data into the univariate_data data table
    univariate_data <- rbind(
        univariate_data,
        covariate_data
    )

}

## merge the univariate and multivariate data
cox_data <- merge(
    univariate_data,
    multivariate_data,
    by = "covariate",
    all.x = T,
    sort = F
)

## save the table into a file in the directory
fwrite(
    cox_data,
    paste(
        snakemake@output[["descriptive_table_dir"]],
        "/",
        "cox_regression",
        ".csv",
        sep = ""
    ),
    sep = ","
)



