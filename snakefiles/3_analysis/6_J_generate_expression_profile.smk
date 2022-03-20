

"""
    Visualize the neurotransmission expression profile for each cluster
"""

rule generate_expression_profile:
    input: 
        genes_highExpressed = "data/others/transformed_expression_data/TCGA_genes_highExpressed.csv",
        kegg_pathway_genes_directory = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes_directory,
        normalized_expression_data_expressionFiltered = "data/expression_data/TCGA_normalized_expression_data_expressionFiltered.csv",
        cluster_group = "data/others/TCGA_IDHall_clusters.csv",
        DE_genes_dir = "data/DEGene_analysis/results_TCGA_IDHall"


    output:
        neurotransmission_expression_profile_dir = directory("data/figures/TCGA_IDHall_neurotransmission_expression_profile")

    params: 
        utils_R = "scripts/r_modules/utils.R",
        pval_cutoff = 0.001,
        logfoldchange_cutoff = 1

    conda:
        "../../envs/r.yaml"

    script: 
        "../../scripts/r_modules/generate_expression_profile.R"