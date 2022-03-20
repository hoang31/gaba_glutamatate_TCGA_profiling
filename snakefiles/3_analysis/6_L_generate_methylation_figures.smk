
"""
    Rule for generating the figures associated with the methylation data
"""

rule generate_methylation_figures:
    input:
        DMgenes_analysis_dir = rules.analyze_methylation_data.output.DMgenes_analysis_dir,
        methylation_information = "data/methylation_data/methylation_information.csv",
        kegg_pathway_genes_directory = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes_directory,
        genes_highExpressed = "data/others/transformed_expression_data/TCGA_genes_highExpressed.csv",


    output:
        methylation_figure_dir = directory("data/figures/TCGA_IDHall_methylation"),
    
    params:
        utils = "scripts/r_modules/utils.R",
        pval_cutoff = 0.001,
        logfoldchange_cutoff = 2

    conda:
        "../../envs/r.yaml"

    script: 
        "../../scripts/r_modules/generate_methylation_figures.R"


