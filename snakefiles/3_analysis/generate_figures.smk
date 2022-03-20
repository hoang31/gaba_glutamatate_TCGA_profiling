
######################################################

"""
            Generate the figures from the epidemiological data
"""

#######################################################

rule generate_figures:
    input:
        epidemiology_data = rules.analyse_epidemiological_data.output.epidemio_dir,
        cluster_group = "data/others/{db}_{data_type}_clusters.csv",
        subgroup_epidemiological_data_all = rules.analyse_epidemiological_data_subgroup.output.subgroup_epidemiological_data_all

    output:
        generated_figures = directory("data/figures/{db}_{data_type}_alteration_figures/")

    conda:
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/generate_figures.R"


rule generate_figures2:
    input:
        subgroup_epidemiological_data_all = "data/TCGA_IDHall_epidemiology_subgroup/all_data.csv"

    output:
        generate_figures2 = "data/generate_figures2"

    conda:
        "../../envs/py_env.yaml"

    script:
        "../../scripts/python_modules/generate_figures.py"
