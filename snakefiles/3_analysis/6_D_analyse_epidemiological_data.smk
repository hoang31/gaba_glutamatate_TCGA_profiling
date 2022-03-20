#########################################################
#########################################################

"""
        ANALYSE EPIDEMIOLOGICAL AND CLINICAL DATA
"""


#########################################################
#########################################################


rule analyse_epidemiological_data:
    input: 
        #subcluster_group = "data/others/{db}_{data_type}_subclusters.csv",
        cluster_group = "data/others/{db}_{data_type}_clusters.csv",
        clinical_data = "data/clinical_data/{db}_clinical_data_cancer",
        external_clinical_data = "data/clinical_data/TCGA_clinical_EXTERNAL.csv",
        purity_score_tumor = "data/others/{db}_purity_score_tumor.csv",
        cnv_data = 'data/others/{db}_{data_type}_alteration_data.csv',

    output: 
        epidemio_dir = directory("data/{db}_{data_type}_epidemiology"),
        epidemiological_data_all = "data/{db}_{data_type}_epidemiology/all_data.csv",
        epidemiological_data_formated = "data/{db}_{data_type}_epidemiology/formated_data.csv",
        epidemiological_data_formated_subset = "data/{db}_{data_type}_epidemiology/epidemiological_data_formated_subset.csv",

    params: 
        wcard = "{data_type}",
        utils_R = "scripts/r_modules/utils.R",

    conda:
        "../../envs/r.yaml"

    script: 
        "../../scripts/r_modules/analyse_epidemiological_data.R"