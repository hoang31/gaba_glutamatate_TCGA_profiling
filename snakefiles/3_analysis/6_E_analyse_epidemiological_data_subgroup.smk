#########################################################
#########################################################

"""
        ANALYSE EPIDEMIOLOGICAL AND CLINICAL DATA of each histological entity of each subgroup
"""


#########################################################
#########################################################


rule analyse_epidemiological_data_subgroup:
    input: 
        epidemiological_data_all = rules.analyse_epidemiological_data.output.epidemiological_data_all

    output: 
       subgroup_epidemiological_data_formated = "data/{db}_{data_type}_epidemiology_subgroup/formated_data.csv",
       subgroup_epidemiological_data_all = "data/{db}_{data_type}_epidemiology_subgroup/all_data.csv"

    params: 
        wcard = "{data_type}",
        utils_R = "scripts/r_modules/utils.R",

    conda:
        "../../envs/r.yaml"

    script: 
        "../../scripts/r_modules/analyse_epidemiological_data_subgroup.R"