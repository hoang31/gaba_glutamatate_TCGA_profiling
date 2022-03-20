
# -*- coding: utf-8 -*-

##########################################################

"""
This script will permit to clean the data from the IGAP data
"""

##########################################################

##### import libraries
import os
from os import listdir
import glob

import pandas as pd
import numpy as np
import re

##########################################################

##### load data

clinical_data_TCGA = pd.read_csv(snakemake.input["clinical_data_TCGA"],  sep = "\t")
id_data_TCGA = pd.read_csv(snakemake.input["id_data_TCGA"],  sep = "\t")

# clinical_data_TCGA = pd.read_csv("/home/hoang/Documents/PROJECTS/gaba_and_glutamate_pathways_in_glioma/data/TCGA/TCGA_clinical_data.tsv",  sep = "\t")
# id_data_TCGA = pd.read_csv("/home/hoang/Documents/PROJECTS/gaba_and_glutamate_pathways_in_glioma/data/TCGA/id_samples.tsv",  sep = "\t")

##########################################################


variables = [
            "File Name",
            "Project ID",
            "Case ID", # corresponding to the the primary key
            "Project ID",
            "Sample ID",
            "Sample Type"
            ]

## extract the specific columns
id_data_TCGA = id_data_TCGA.loc[:,variables]

## replace the blank with "_"
variables = [variable.replace(" ", "_")
            for variable in variables
            ]

## rename the column name
id_data_TCGA.columns = variables


##########################################################


##### Clinical data

variables = [
            "submitter_id", # corresponding to the primary key
            "age_at_index",
            "days_to_death",
            "gender",
            "vital_status",
            #"age_at_diagnosis",
            "primary_diagnosis",
            "site_of_resection_or_biopsy",
            "year_of_diagnosis",
            ]

## extract the data associated with the specific variables
clinical_data_TCGA = clinical_data_TCGA.loc[:,variables]


TCGA_clinical_data = pd.merge(id_data_TCGA, clinical_data_TCGA,
                how = "left",
                left_on = "Case_ID",
                right_on = "submitter_id",
                )

TCGA_clinical_data = TCGA_clinical_data.drop_duplicates()


# TCGA_clinical_data.to_csv(path_or_buf = "/home/hoang/Desktop/TCGA_test2.tsv", sep = "\t", index = None, header=True)

##########################################################

##### Removing the replicates from the data

## retreive the cases having replicates
cases_with_replicates = pd.DataFrame(TCGA_clinical_data["Case_ID"].value_counts())
cases_with_replicates = cases_with_replicates[cases_with_replicates["Case_ID"] > 2].index

## retreinve the ID samples from the cases having replicates
samples_with_replicates =list( TCGA_clinical_data.loc[TCGA_clinical_data["Case_ID"].isin(cases_with_replicates)]["Sample_ID"])

## replicate to keep
replicates_to_kept = [] ## initialize the list that will contain all the replicates to keep
for i in range(0,len(cases_with_replicates)) :

    ## extract the replicate id with the case id
    filt = list(filter(lambda x: re.search(cases_with_replicates[i], x),  samples_with_replicates))

    ## extract the replicate having a "A" in the id replicate
    filt2 = list(filter(lambda x: re.search(r'[0-5][0-9]A', x),  filt))

    ## take the higher replicate id number
    higher_index = max([part[13:15] for part in filt2])

    ## keep the replicate with a A and the higher index
    replicat_to_kept = cases_with_replicates[i] + "-" + higher_index + "A"
    replicates_to_kept.append(replicat_to_kept)

## remove case with replicates
TCGA_clinical_data_final = TCGA_clinical_data[TCGA_clinical_data["Case_ID"].isin(cases_with_replicates) == False]
TCGA_clinical_data_final.append(TCGA_clinical_data[TCGA_clinical_data["Sample_ID"].isin(replicates_to_kept)])

## In this TCGA data, we have the case "TCGA-06-0156" having two replicates but with the same ID samples. So we have to remove manually keeping one of two.
## We choose to keep the data from the : 18fbf794-a85b-4b0d-9a5c-c1c0e0ede3d8.htseq.counts.gz and 18fbf794-a85b-4b0d-9a5c-c1c0e0ede3d8.FPKM.txt.gz. We have to remove the c341aa36-a431-477a-9fae-4550c4eea047.FPKM.txt.gz and c341aa36-a431-477a-9fae-4550c4eea047.htseq.counts.gz files.

files_to_remove = ["c341aa36-a431-477a-9fae-4550c4eea047.FPKM.txt.gz",
                "c341aa36-a431-477a-9fae-4550c4eea047.htseq.counts.gz"]

TCGA_clinical_data_final = TCGA_clinical_data_final.loc[TCGA_clinical_data_final["File_Name"].isin(files_to_remove) == False]



##########################################################

## add case associated with the replicates to keep
TCGA_clinical_data_final = TCGA_clinical_data_final.append(TCGA_clinical_data[TCGA_clinical_data["Sample_ID"].isin(replicates_to_kept)])

## In this TCGA data, we have the case "TCGA-06-0156" having two replicates but with the same ID samples. So we have to remove manually keeping one of two.
## We choose to keep the data from the : 18fbf794-a85b-4b0d-9a5c-c1c0e0ede3d8.htseq.counts.gz and 18fbf794-a85b-4b0d-9a5c-c1c0e0ede3d8.FPKM.txt.gz. We have to remove the c341aa36-a431-477a-9fae-4550c4eea047.FPKM.txt.gz and c341aa36-a431-477a-9fae-4550c4eea047.htseq.counts.gz files.

files_to_remove = ["c341aa36-a431-477a-9fae-4550c4eea047.FPKM.txt.gz",
                "c341aa36-a431-477a-9fae-4550c4eea047.htseq.counts.gz"]

TCGA_clinical_data_final = TCGA_clinical_data_final.loc[TCGA_clinical_data_final["File_Name"].isin(files_to_remove) == False]


##########################################################


##### Healthy and cancer samples

TCGA_healthy_cases = TCGA_clinical_data_final[TCGA_clinical_data_final["Sample_Type"] == "Solid Tissue Normal"]

TCGA_cancer_cases = TCGA_clinical_data_final[TCGA_clinical_data_final["Sample_Type"] != "Solid Tissue Normal"]


print(TCGA_healthy_cases.shape)
print(TCGA_cancer_cases.shape)

# print(TCGA_healthy_cases.shape)
# print(TCGA_cancer_cases.shape)
# print(TCGA_clinical_data_final.shape)


##########################################################


##### Expression data
TCGA_expression_counts = snakemake.input["expression_data_TCGA_count"]
TCGA_expression_fpkm = snakemake.input["expression_data_TCGA_fpkm"]

# TCGA_expression_counts = "/home/hoang/Documents/PROJECTS/gaba_and_glutamate_pathways_in_glioma/data/TCGA/expression_count.tsv"
# TCGA_expression_FPKM = "/home/hoang/Documents/PROJECTS/gaba_and_glutamate_pathways_in_glioma/data/TCGA/expression_FPKM.tsv"

##### function for cleaning the data
## expressiom data correspond to the the expression data with in column each sample and in row the gene IDs
## the expression type correspond to "count" or "FPKM"
## case_type corresponding to : "normal" or "cancer"

def clean_data_expression(expression_data_path, expression_type, case_type) :

    ## load csv data
    TCGA_expression_data = pd.read_csv(expression_data_path, sep = "\t")

    ## remove the column that are replicated
    TCGA_expression_data = TCGA_expression_data[TCGA_expression_data.columns.drop(list(TCGA_expression_data.filter(regex='genes.[1-9][0-9]{0,2}')))]


    if case_type == "cancer" :

        ## extract the file names from the id sample data
        file_names = ['.'.join(file_name.split(".")[0:-1])
                    for file_name in TCGA_cancer_cases["File_Name"].tolist()]

        ## retreive the name file of FPKM and COUNT file, adding the variable name "genes"
        files = list(filter(lambda x: re.search(expression_type, x),  file_names))
        files = ["genes"] + files

    else :

        ## extract the file names from the id sample data
        file_names = ['.'.join(file_name.split(".")[0:-1])
                    for file_name in TCGA_healthy_cases["File_Name"].tolist()]

        ## retreive the name file of FPKM and COUNT file, adding the variable name "genes"
        files = list(filter(lambda x: re.search(expression_type, x),  file_names))
        files = ["genes"] + files



    ## add the extension .gz for extract from the clinical data, the case ID for putting as the column name
    file_gz = [filename + ".gz" for filename in files]
    # case_id = TCGA_clinical_data_final.loc[TCGA_clinical_data_final["File_Name"].isin(file_gz)].loc[:,["Case_ID"]]
    case_id = TCGA_clinical_data_final.loc[TCGA_clinical_data_final["File_Name"].isin(file_gz)]["Case_ID"].tolist()

    print(type(case_id))
    case_id = ["genes"] + case_id


    ## extract the expression data associated with our files
    TCGA_expression_data = TCGA_expression_data.loc[:,files]

    TCGA_expression_data.columns = case_id

    return(TCGA_expression_data)


##########################################################


##### Write the data

## clincial data
TCGA_cancer_cases.to_csv(path_or_buf = snakemake.output["TCGA_clinical_data_cancer"], sep = ",", index = None, header=True)
TCGA_healthy_cases.to_csv(path_or_buf = snakemake.output["TCGA_clinical_data_normal"], sep = ",", index = None, header=True)

## expression data
print("FPKM CANCER")
clean_data_expression(TCGA_expression_fpkm, "FPKM", "cancer").to_csv(path_or_buf = snakemake.output["TCGA_expression_data_fpkm_cancer"], sep = ",", index = None, header=True)
print("FPKM NORMAL")
clean_data_expression(TCGA_expression_fpkm, "FPKM", "normal").to_csv(path_or_buf = snakemake.output["TCGA_expression_data_fpkm_normal"], sep = ",", index = None, header=True)

print("COUNT CANCER")
clean_data_expression(TCGA_expression_counts, "count", "cancer").to_csv(path_or_buf = snakemake.output["TCGA_expression_data_count_cancer"], sep = ",", index = None, header=True)
print("COUNT NORMAL")
clean_data_expression(TCGA_expression_counts, "count", "normal").to_csv(path_or_buf = snakemake.output["TCGA_expression_data_count_normal"], sep = ",", index = None, header=True)



##########################################################
