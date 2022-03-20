
#########################################################
#########################################################

"""
        SURVIVAL ANALYSIS of each clusters
"""

#########################################################
#########################################################


rule calculate_survival:
    input: 
        cluster_group = "data/others/{db}_{data_type}_clusters.csv",
        subgroup_epidemiological_data_all = 'data/TCGA_IDHall_epidemiology_subgroup/all_data.csv',
        clinical_data = "data/clinical_data/{db}_clinical_data_cancer",
        clinical_data_full = "data/clinical_data/TCGA_clinical_FULL.tsv",
        TCGA_clinical_EXTERNAL = "data/clinical_data/TCGA_clinical_EXTERNAL.csv",
        epidemiological_data_all = "data/TCGA_IDHall_epidemiology_subgroup/all_data.csv"

    output: 
        survival_curve_dir = directory("data/figures/{db}_{data_type}_survival_curves")

    params: 
        wcard = "{data_type}",
        utils_R = "scripts/r_modules/utils.R",

    conda:
        "../../envs/r.yaml"

    script: 
        "../../scripts/r_modules/calculate_survival.R"