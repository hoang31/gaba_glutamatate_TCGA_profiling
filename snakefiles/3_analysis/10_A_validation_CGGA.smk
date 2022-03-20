
######################################################

"""
        Validation with CGGA database
"""

#######################################################


rule validation_CGGA:
    input:
        expression_data_CGGA = 'data/expression_data/CGGA_expression_data.csv',
        clinical_data_CGGA = 'data/clinical_data/CGGA_clinical_data.csv',
        gene_of_interest = 'data/others/transformed_expression_data/TCGA_genes_highExpressed.csv',
        cgga_immune_estimation_dir = 'data/CGGA_timer2_estimation',
        kegg_pathway_genes_directory = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes_directory,


    output:
        cgga_validation_dir = directory("data/figures/CGGA_validation")

    params:
        heatmap_function = "scripts/r_modules/heatmap_generation.R",
        utils = "scripts/r_modules/utils.R",
        clustering_metrics = "scripts/r_modules/clustering_metrics.R",


    conda:
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/validation_CGGA.R"
