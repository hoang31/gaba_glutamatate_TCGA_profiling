
# -*- coding: utf-8 -*-

##########################################################

"""
This script will permit to clean the data from the IGAP data
"""

##########################################################

##### import libraries
import os
import pandas as pd
import numpy as np
import re

##########################################################

##### load data
expression_data_IGAP = pd.read_csv(snakemake.input["expression_data_IGAP"], sep = "\t", skipinitialspace = True)

clinical_data_IGAP = pd.read_csv(snakemake.input["clinical_data_IGAP"],  sep = "\t", skipinitialspace = True)
# clinical_data_IGAP.fillna("TO_DELETE", inplace=True) ## fill the NA values for deleting

print(expression_data_IGAP.columns)

##########################################################


##### clinical data

## variable to select
variables = ["Allen Institute #",
            "Histopathology", ## = histology + gradde
            "Gender",
            "KPS",
            "Age At Initial Diagnosis (in years)",
            "Survival (in days)",
            "Cause Of Death",
            "Surgery (primary or recurrent)",

            ## biomarkers
            "1p19q_deletion",
            "EGFR",
            "PTEN",
            "MGMT PCR",
            "MGMT",
            "IDH1",
            ]

## exctract the variable data from the clinical data and create the merged clinical data
IGAP_clinical_data = clinical_data_IGAP.loc[:,variables]

## rename columns
colnames = [colname.replace(" ", "_")
            for colname in IGAP_clinical_data.columns]

IGAP_clinical_data.columns = colnames


##########################################################


##### expression data

## extract the samples from the expression data having clinical data
id_samples = expression_data_IGAP.columns.intersection(list(IGAP_clinical_data["Allen_Institute_#"])).tolist()

## add the gene symbols column
id_samples2 = ["genesymbol"] + id_samples

## extract the expression related to these id_samples
expression_data_IGAP = expression_data_IGAP.loc[:,id_samples2]

## change the column name
expression_data_IGAP.columns = ["Gene_Name"] + id_samples

##### clinical data

## extrac the clinical data related to these id_samples
IGAP_clinical_data = IGAP_clinical_data[IGAP_clinical_data["Allen_Institute_#"].isin(id_samples2)]



##########################################################


##### change values for standardization with other databases

## OS
IGAP_clinical_data["Cause_Of_Death"] =  IGAP_clinical_data["Cause_Of_Death"].replace("Deceased due to tumor progression", "dead")
IGAP_clinical_data["Cause_Of_Death"].fillna("alive", inplace=True)

## histology and grade
histology = []
grade = []
for i in IGAP_clinical_data["Histopathology"] :
    match = re.search("Glioblastoma", i)
    if match :
        histology.append("GBM")
        grade.append("WHO IV")
    else :
        histology.append(i)
        grade.append("Unknown")

IGAP_clinical_data["Histology"] = histology
IGAP_clinical_data["Grade"] = grade
IGAP_clinical_data = IGAP_clinical_data.drop(columns = ["Histopathology"])


## rename the variables
IGAP_clinical_data = IGAP_clinical_data.rename(columns={
                                                        "Allen_Institute_#" : "sample_id",
                                                        "Cause_Of_Death": "vital_status",
                                                        "Age_At_Initial_Diagnosis_(in_years)" : "Age_in_years",
                                                        "Survival_(in_days)" : "OS",
                                                        "Surgery_(primary_or_recurrent)" : "reccurence",
                                                        "1p19q_deletion" : "1p19q_deletion_status",
                                                        "IDH1" : "IDH_mutation_status",
                                                        "EGFR" : "EGFR_status",
                                                        "PTEN" : "PTEN_status"
                                                        }
                                                )

## add the database variable
IGAP_clinical_data["Database"] = (["IGAP" for i in range(0, IGAP_clinical_data.shape[0])])


##########################################################

##### write the data : clinical and expression data

expression_data_IGAP.to_csv(path_or_buf = snakemake.output["IGAP_expression_data"], sep = ",", index = None, header=True)
IGAP_clinical_data.to_csv(path_or_buf = snakemake.output["IGAP_clinical_data"], sep = ",",index = None, header=True)
