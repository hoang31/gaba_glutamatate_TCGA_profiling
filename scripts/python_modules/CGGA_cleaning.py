
# -*- coding: utf-8 -*-

##########################################################

"""
This script will permit to clean the data from the CGGA data
"""

##########################################################

##### import libraries
import os
import pandas as pd
import numpy as np
import re

##########################################################

##### load data
expression_data_CGGA1 = pd.read_csv(snakemake.input["expression_data_CGGA1"], sep = "\t")
expression_data_CGGA2 = pd.read_csv(snakemake.input["expression_data_CGGA2"],  sep = "\t")
clinical_data_CGGA1 = pd.read_csv(snakemake.input["clinical_data_CGGA1"], sep = "\t|\s\t|\s\s\t", engine = "python")
clinical_data_CGGA2 = pd.read_csv(snakemake.input["clinical_data_CGGA2"], sep = "\t|\s\t|\s\s\t", engine = "python")


##########################################################
##### clinical data

## variable to select
variables = ["CGGA_ID",
            "Histology",
            "Grade",
            "Gender",
            "Age",
            "OS",
            "Censor",
            "IDH_mutation_status",
            "1p19q_codeletion_status"
            ]

## exctract the variable data from the clinical data and create the merged clinical data
CGGA_clinical_data = clinical_data_CGGA1.loc[:,variables].append(clinical_data_CGGA2.loc[:,variables]).reset_index(drop=True)
CGGA_clinical_data.fillna("Unknown", inplace=True) ## remplace the NA values

## rename column names
CGGA_clinical_data = CGGA_clinical_data.rename(columns = {
                                                        "CGGA_ID" : "sample_id",
                                                        "Age" : "Age_in_years",
                                                        "Censor" : "vital_status"
                                                        }
                                                )
## change the variable
CGGA_clinical_data["vital_status"] =  CGGA_clinical_data["vital_status"].replace(0, "alive")
CGGA_clinical_data["vital_status"] =  CGGA_clinical_data["vital_status"].replace(1, "dead")

## add the database variable
CGGA_clinical_data["Database"] = (["CGGA" for i in range(0, CGGA_clinical_data.shape[0])])

## extract the recurrence of the tumors from the histology column
reccurence = []
histology = []
for i in range(0,(CGGA_clinical_data.shape[0])) :
    # print("................................")
    # print(i)
    # print(CGGA_clinical_data.loc[CGGA_clinical_data.index[i], "Histology"])
    match = re.split(pattern = "r", string = CGGA_clinical_data.loc[CGGA_clinical_data.index[i], "Histology"])

    if (len(match) == 1) :
        if match[0] == "Unknown" :
            reccurence.append("Unknown")
            histology.append("Unknown")
        else :
            reccurence.append("Primary")
            histology.append(match[0])
    else :
        reccurence.append("Recurrent")
        histology.append(match[1])

CGGA_clinical_data.drop(columns = ["Histology"]) ## remove the old Histology column

## adding new colunm in the data
CGGA_clinical_data["Histology"] = histology
CGGA_clinical_data["reccurence"] = reccurence


##########################################################


##### expression data

CGGA_expression_data = pd.merge(expression_data_CGGA1, expression_data_CGGA2,
        how = "inner",
        on = "Gene_Name",
        )

##########################################################

##### write the data : clinical and expression data

CGGA_expression_data.to_csv(path_or_buf = snakemake.output["CGGA_expression_data"], sep = ",",index = None, header=True)
CGGA_clinical_data.to_csv(path_or_buf = snakemake.output["CGGA_clinical_data"], sep = ",", index = None, header=True)
