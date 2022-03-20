
"""
    Rule for generate the figure of correlation between methylation data and gene expression data
"""

rule generate_correlation_methylation_expression_figures:
    input:
        DMgenes_analysis_dir = rules.analyze_methylation_data.output.DMgenes_analysis_dir,
        methylation_information = "data/methylation_data/methylation_information.csv",
        kegg_pathway_genes_directory = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes_directory,
        genes_highExpressed = "data/others/transformed_expression_data/TCGA_genes_highExpressed.csv",
        expression_data = "data/expression_data/TCGA_normalized_expression_data_expressionFiltered.csv"


    output:
        methylation_figure_correlation_dir = directory("data/figures/TCGA_IDHall_methylation_correlation"),
        methylation_expression_data = "data/others/gene_of_interest_data/1_methylation_expression_data.csv"
    
    params:
        utils = "scripts/r_modules/utils.R",
        #pval_cutoff = 0.001,
        #logfoldchange_cutoff = 2

    conda:
        "../../envs/r.yaml"

    script: 
        "../../scripts/r_modules/generate_correlation_methylation_expression_figures.R"

