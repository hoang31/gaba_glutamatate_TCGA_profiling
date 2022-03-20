
"""
        Analyse the methylation data
"""

rule analyze_methylation_data:
    input:
        methylation_information = "data/methylation_data/methylation_information.csv",
        methylation_data = "data/methylation_data/methylation_data.csv",
        genes_highExpressed = "data/others/transformed_expression_data/TCGA_genes_highExpressed.csv",
        cluster_data = "data/others/TCGA_IDHall_clusters.csv",
        kegg_pathway_genes_directory = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes_directory,

    output:
        DMgenes_analysis_dir = directory("data/DMgenes_analysis"),
    
    params:
        utils = "scripts/r_modules/utils.R"

    conda:
        "../../envs/r.yaml"

    script: 
        "../../scripts/r_modules/analyze_methylation_data.R"


