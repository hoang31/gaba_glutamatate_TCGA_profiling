
"""
    Generation of the tables that describe the data
"""

rule generate_tables:
    input:
        formated_epidemiological_data = "data/TCGA_IDHall_epidemiology/epidemiological_data_formated_subset.csv",
        survival_data_dir = "data/figures/TCGA_IDHall_survival_curves"
    
    output:
        descriptive_table_dir = directory("data/figures/TCGA_IDHall_descriptive_tables")
    
    conda:
        "../../envs/r.yaml"

    script: 
        "../../scripts/r_modules/generate_tables.R"