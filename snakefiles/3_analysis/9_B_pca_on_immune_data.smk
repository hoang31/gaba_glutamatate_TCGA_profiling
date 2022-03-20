
"""
        Infering immune composition for each glioma samples using timer2 
"""

rule pca_on_immune_data:
    input:
        immune_composition_data = rules.download_immune_data_TCGA.output.TCGA_immune_data,
        clinical_data = 'data/clinical_data/TCGA_clinical_data_cancer',
        cluster_group = "data/others/TCGA_IDHall_clusters.csv",

    output:
        #pca_figure = 'data/pca.svg'
        pca_figure_dir = directory("data/figures/TCGA_IDHall_immune_pca")

    params:
        utils_R = "scripts/r_modules/utils.R"
    conda:
        "../../envs/r_pca.yaml"
    script:
        "../../scripts/r_modules/pca_on_immune_data.R" 