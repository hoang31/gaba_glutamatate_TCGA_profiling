
"""
        Infering immune composition for each glioma samples using timer2 
"""

rule infer_immune_type_composition:
    input:
        immune_composition_data = rules.download_immune_data_TCGA.output.TCGA_immune_data,
        clinical_data = 'data/clinical_data/TCGA_clinical_data_cancer',
        cluster_group = "data/others/TCGA_IDHall_clusters.csv",

    output:
        immune_composition_figures_dir = directory('data/figures/TCGA_IDHall_immune_composition_figures')


    params:
        utils_R = "scripts/r_modules/utils.R"
    conda:
        "../../envs/r.yaml"
    script:
        "../../scripts/r_modules/infer_immune_type_composition.R" 