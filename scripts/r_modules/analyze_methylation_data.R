
#####################################
## Analyze the methylation data associated with the gene of interest
#####################################


#####################################
## Load the libraries
#####################################


library(data.table)
library(tidyverse)
library(RColorBrewer)


#####################################
## load the data 
#####################################


## read the data table that contain the methylation information
methylation_information <- fread(
    snakemake@input[["methylation_information"]],
    sep = ","
)

## read the data table that contain the methylation data
methylation_data <- fread(
    snakemake@input[["methylation_data"]],
    #drop = c(3:680),
    sep = ","
)

#methylation_data <- fread(
#    snakemake@input[["methylation_data"]],
#    #drop = c(3:680),
#    sep = ","
#)[1:50,]

## load the genes related to GABA glutamate and calcium pathway
gene_vector <- unique(
    fread(
        snakemake@input[["genes_highExpressed"]],
        sep = ","
    )[, gene_name]
)

## load the cluster data
cluster_data <- fread(
    snakemake@input[["cluster_data"]],
    sep = ","

)

## add the sample healthy into the cluster_data
healthy_sample_vector <- c(
    "TCGA-06-AABW",
    "TCGA-74-6573"
)
healthy_cluster_data <- data.table(
    "sample_id" = healthy_sample_vector,
    "cluster" = 'healthy'
)

## rbind with the cluster data
cluster_data <- rbind(
    cluster_data,
    healthy_cluster_data
)

## load the utils functions
source(snakemake@params[["utils"]])

## color used for the figures
color_palette <- brewer.pal(9, "Set1")

## create the directory that will contain the data
dir.create(snakemake@output[["DMgenes_analysis_dir"]])


#####################################
## Cleaning the methylation data (step that is not necessary)
#####################################

##########
## Merge the methylation data and write it into files, by bins
##########


### create a directory that will contain the bin methylation data
#dir.create("data/methylation_data/all")

### retreive the file path of all the methylation data
#methylation_path_vector <- list.files(
#    "data/methylation_data",
#    pattern = "Methylation450",
#    full.names = T,
#    recursive = T
#)

### create the bins
#bins_vector <- c(seq(1, length(methylation_path_vector), 100),length(methylation_path_vector))

## for all the bins, retreive methylation data and merge it into a file
##for (i_bins in seq(1, (length(bins_vector) - 1))) {
##for (i_bins in seq(1, 5)) {
#for (i_bins in seq(6, (length(bins_vector) - 1))) {

#    ## generate the inferior and superior borns
#    bins_inf <- bins_vector[i_bins]
#    bins_sup <- (bins_vector[i_bins+1]) - 1

#    ## retreive the methylation data name related to the bins
#    methylation_path_vector_subset <- methylation_path_vector[bins_inf:bins_sup]

#    ## for each methylation file, read it, clean it
#    for (i_methylation_data in seq(1, length(methylation_path_vector_subset))) {
        
#        print(i_methylation_data)
#        ## extract the path of one file
#        file_path <- methylation_path_vector_subset[i_methylation_data]
        
#        ## retreive the id of the sample
#        sample_id <- str_match(file_path, "lvl-3.(.*?).gdc_hg38.txt")[2]
#        sample_id <- substr(sample_id, start = 1, stop = 20)

#        ## read the file, extract the column of interest and change the column name by the sample id
#        methylation_data <- fread(
#            file_path,
#            sep = "\t"
#        )[, c("Composite Element REF", "Beta_value")][, setnames(.SD, "Beta_value", sample_id)]

#        ## if we are for the first methylation data, initialize the data tble that contain the merged data with it
#        if (i_methylation_data == 1) {
#            methylation_data_dt <- copy(methylation_data)
#            next
#        }

#        ## merge with all the data together
#        methylation_data_dt <- merge(
#            methylation_data_dt,
#            methylation_data,
#            by = "Composite Element REF",
#            sort = F
#        )

#        ## remove the na values
#        #methylation_data_dt <- na.omit(methylation_data_dt)

#    }

#    ## generate the path of the file
#    file_path <- paste(
#        "data/methylation_data/all",
#        bins_inf,
#        "_",
#        bins_sup,
#        ".csv",
#        sep = ""
#    )

#    ## write the merged methylation data
#    fwrite(
#        methylation_data_dt,
#        file = file_path,
#        sep = ","
#    )
#}

#break

##########
## Merge each bin data together into a same files
##########

### load the name files 
#files <- list.files(
#    "data/methylation_data/all",
#    pattern = ".csv",
#    full.names = T
#)

### read the first file to initialize the data table that will contain the merged methylation data
#merge_data <- fread(files[1], sep = ",")

### for all the other files, merge it with the merge_data data table
#for (i_file in seq(2, length(files))) {
    
#    ## retreive file path
#    file_path <- files[i_file]

#    ## read file using the file path
#    file_dt <- fread(file_path, sep = ",")

#    ## merge the data into the merge_data data table
#    merge_data <- merge(
#        merge_data,
#        file_dt,
#        by = "Composite Element REF",
#        all.x = T,
#        sort = F
#    )
#}

### write the data into file
#fwrite(merge_data, "data/methylation_data/methylation_data.csv")

#break

##########
## Create the table that will contain the genes related to each probe
##########

## create the table that will contain each propes associated with the related genes
#methylation_information <- fread(
#    "data/methylation_data/GBM/0000c40e-9d45-4446-9dd9-a4676224d0ce/jhu-usc.edu_GBM.HumanMethylation450.7.lvl-3.TCGA-19-5955-01A-11D-1697-05.gdc_hg38.txt",
#    select = c("Composite Element REF", "Gene_Symbol"),
#    sep = "\t"
#)

### extract the column of interest from the methylation information
#methylation_information <- methylation_information[, c("Composite Element REF", "Gene_Symbol"), with = F]

### extract the probe name into a vector
#probe_name_vector <- unlist(methylation_information[, "Composite Element REF", with = F])

### initialize the table that will contain the probe and gene informaiton
#methylation_information_interest <- data.table(
#    "probe_id" = character(),
#    "gene_name" = character(),
#    "probe_name" = character()
#)

#print(length(probe_name_vector))
### for each probe, retreive the names of the genes that are related to the probes
#for (i_probe in seq(1, length(probe_name_vector))) {

#    print(i_probe)
#    ## identify the probe names 
#    probe_name <- probe_name_vector[i_probe]

#    ## identify the unique genes associated with the probe name
#    gene_name_vector <- unique(
#            unlist(
#            str_split(
#                string = unlist(methylation_information[`Composite Element REF` == probe_name, Gene_Symbol]),
#                pattern = ";"
#            )
#        )
#    )

#    ## generate the table that contain the probe name and the gene names
#    subset_dt <- data.table(
#        "probe_id" = probe_name,
#        "gene_name" = gene_name_vector
#    )

#    ## add the probe name
#    subset_dt[, probe_name := paste(probe_id, gene_name, sep = "_")]

#    ## merge the subset_dt with the mutlation_data
#    methylation_information_interest <- rbind(
#        methylation_information_interest,
#        subset_dt
#    )

#}

### remove the rows that does not contain gene name associated with the probes
#methylation_information_interest <- methylation_information_interest[gene_name != ".", ]

### write the data into csv file
#fwrite(
#    methylation_information_interest,
#    "data/all/methylation_information.csv",
#    sep = ","
#)


#####################################
## Analysis
#####################################


##########
## extract the methylation data of interest
##########

## modification of the methylation information
methylation_information[grep(gene_name, pattern = 'SLC25A6'), gene_name := substr(gene_name, start = 1, stop = 7)]

## extract the methylation probe id associated with the gene of interest
methylation_information <- methylation_information[gene_name %in% gene_vector, ][order(gene_name),]

## extract the probe names associated with the gene of interest
probe_name_of_interest_vector <- unlist(methylation_information[, probe_id])

## extract the methylation data associated with the probe names of interest
methylation_data <- methylation_data[`Composite Element REF` %in% probe_name_of_interest_vector,]

## rename the column that contain the probes id
methylation_data <- methylation_data[, setnames(.SD, "Composite Element REF", "probe_id")]


##########
## clean the methylation data
##########


## extract the methylation information of the probes of interest
methylation_information_of_interest <- methylation_information[probe_id %in% probe_name_of_interest_vector,]

## extract the methylation information from the probes of interest
methylation_information_counting <- methylation_information_of_interest[, .(counting = .N), by = "probe_id"][order(counting),]

print(dim(methylation_data))
print(nrow(methylation_information_of_interest))
print(nrow(methylation_information_counting))


## extract the probes that have different names
probes_vector <- unlist(methylation_information_counting[counting > 1, probe_id])


## initialize the data table that will contain the methylation data of the probes that are associated with several genes
duplicated_methylation_data <- copy(methylation_data[0,])

## for each probes, we will extract the methylation data and put the new id
for (i_probe in seq(1, length(probes_vector))) {

    ## retreive the name of the probe
    probe_of_interest_name <- probes_vector[i_probe]
    #print(probe_of_interest_name)

    ## retreive the information concerning the probe
    meth_information <- methylation_information[probe_id == probe_of_interest_name, ]

    ## retreive the methylation data concerning the probe
    meth_data <- methylation_data[probe_id %in% probe_of_interest_name,]

    ## retreive the number of time that the probe appears
    time_counting <- nrow(meth_information)

    ## dupplicate the meth data x time depending to the time_counting
    meth_data <- meth_data[rep(1, time_counting),]

    ## retreive the probe name for each row
    probe_name_vector <- unlist(meth_information[, probe_name])

    ## rename the rows of the meth_data
    meth_data[, probe_id := probe_name_vector]

    ## put the meth_data to the duplicated methylation data table
    duplicated_methylation_data <- rbind(
        duplicated_methylation_data,
        meth_data
    )

}

## remove the probes associating several genes from the methylation data
methylation_data <- methylation_data[!(probe_id %in% probes_vector),]

## function to transform the probe id to probe name
probe_id_to_probe_name <- function(
    probe_id_input,
    methylation_information_input
) {

    ## retreive the probe name from the probe id
    probe_name_output <- methylation_information_input[probe_id %in% probe_id_input,]

    ## reorder the colum based on the input
    probe_name_output <- unlist(probe_name_output[match(probe_id_input, probe_id), probe_name])


    return(probe_name_output)
}

## rename probe id to probe name
methylation_data[, probe_id := probe_id_to_probe_name(probe_id,methylation_information)]

## add the duplicated methylation to the methylation data
methylation_data <- rbind(
    methylation_data,
    duplicated_methylation_data
)

## do the transposition of the methylation data
methylation_data <- transpose_datatable(
    methylation_data,
    column_name = "probe_id",
    new_name = "sample_id"
)[order(sample_id),]

## split the id 
methylation_data[, sample_id2 := substr(sample_id, start = 1, stop = 12) ]

## remove one of the replicated if there are
duplicated_id_vector <- methylation_data[,sample_id2]
duplicated_id_vector <- unique(duplicated_id_vector[duplicated(duplicated_id_vector)])

## for each deuplicated id, keep only just one replicated
for (i_duplicated_id in seq(1, length(duplicated_id_vector))) {

    ## extract the i-th duplicated id 
    duplicated_id <- duplicated_id_vector[i_duplicated_id]

    ## extract the sample ids associated with this duplicated id
    sample_id_vector <- unlist(methylation_data[sample_id2 == duplicated_id, sample_id])

    ## select the sampples to remove
    sample_id_vector_to_remove <- sample_id_vector[2:(length(sample_id_vector))]

    ## remove the duplicated from the methylation
    methylation_data <- methylation_data[!(sample_id %in% sample_id_vector_to_remove),]

}

## remove one of sample id column and rename the new sampld id column
methylation_data <- methylation_data[, sample_id := NULL][, setnames(.SD, "sample_id2", "sample_id")]

## merge the methylation data with the cluster data
methylation_data <- merge(
    cluster_data,
    methylation_data,
    by = "sample_id"
    #all.x = TRUE,
)

## retreive the probes of interest from the methylation data
probe_name_of_interest_vector <- colnames(methylation_data)[!(colnames(methylation_data) %in% c("sample_id", "cluster"))]



##########
## Do the statistical ananlysis and retreive the site that are differentially methylated between the clusters
##########


## do all the paired combination of cluster
cluster_combination <- combn(
    unique(unlist(cluster_data[,cluster])),
    m = 2
)

## initialize the list that will contain all the combination statistical results
statistical_results_list <- list()

## initialize the list that will contain all significant methylated sites
statistical_probes_list <- list()

## for all the cluster combination
for (i_cluster_combination in seq(1, ncol(cluster_combination))) {
#for (i_cluster_combination in seq(1, 1)) {

    ## initialize the data table that will contain the pvalues of the combination
    combination_statistical_results <- data.table(
        probe_id = character(),
        pvalue = numeric(),
        logfoldchange = numeric()
    )

    ## retreive the combination
    name_condition1 <- cluster_combination[1, i_cluster_combination]
    name_condition2 <- cluster_combination[2, i_cluster_combination]

    ## for each probes, do the statistical analysis
    for (i_probe in seq(1, length(probe_name_of_interest_vector))) {
    #for (i_probe in seq(1, 10)) {

        print(i_probe)

        ## retreive the probe id
        probe_id <- probe_name_of_interest_vector[i_probe]

        ## extract the beta values associated with the probe from the methylation data
        values_condition1 <- as.numeric(unlist(methylation_data[cluster == name_condition1, probe_id, with = F]))
        values_condition2 <- as.numeric(unlist(methylation_data[cluster == name_condition2, probe_id, with = F]))

        if ((sum(is.na(values_condition1)) > (length(values_condition1) - 1)) | (sum(is.na(values_condition2)) > (length(values_condition2) - 1))) {
            
            ## create a data table with the probe id and the pvalue
            subset_dt <- data.table(
                probe_id = probe_id,
                pvalue = NA,
                logfoldchange = NA

            )

            ## merge the data with the combination_statistical_results data table
            combination_statistical_results <- rbind(
                combination_statistical_results,
                subset_dt
            )

            next
        }

        ## calculate the log fold change
        logfoldchange <- log2(mean(values_condition1)/mean(values_condition2))

        ## do the statiscal analysis
        p_val <- wilcox.test(
            values_condition1,
            values_condition2,
        )$p.value

        ## create a data table with the probe id and the pvalue
        subset_dt <- data.table(
            probe_id = probe_id,
            pvalue = p_val,
            logfoldchange = logfoldchange
        )

        ## merge the data with the combination_statistical_results data table
        combination_statistical_results <- rbind(
            combination_statistical_results,
            subset_dt
        )
    }

    ## adjust the pvalue
    combination_statistical_results <- combination_statistical_results[, adj_pvalue := p.adjust(pvalue, method = "bonferroni", n = nrow(combination_statistical_results))][order(adj_pvalue, decreasing = F),]

    ## name of the combnation analysis
    statistical_name <- paste(
        name_condition1,
        "_vs_",
        name_condition2,
        sep = ""
    )

    ## put the results into the list 
    statistical_results_list[[statistical_name]] <- combination_statistical_results

}


######################################
### Write data
######################################


## write the differential methylated gene results
for (i_result in seq(1, length(statistical_results_list))) {

    ## name of the results
    result_name <- names(statistical_results_list[i_result])

    ## extract the data table that contain the results
    data_results <- as.data.table(statistical_results_list[[result_name]])

    ## path of the data
    file_path <- paste(
        snakemake@output[["DMgenes_analysis_dir"]],
        "/",
        result_name,
        ".csv",
        sep = ""
    )

    ## write the data
    fwrite(
        data_results,
        file_path,
        sep = ","
    )
}

## write the data table that contain the methylation data
fwrite(
    methylation_data,
    paste(
        snakemake@output[["DMgenes_analysis_dir"]],
        "/",
        'methylation_data_filtered',
        ".csv",
        sep = ""
    )
)
