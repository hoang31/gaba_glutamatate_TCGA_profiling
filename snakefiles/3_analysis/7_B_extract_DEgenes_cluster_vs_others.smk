
######################################################

"""
            Analyse the differential expressed genes
"""

#######################################################


rule extract_DEgenes_cluster_vs_others:
    input:
        DEgene_results_dir = "data/DEGene_analysis_cluster_vs_others",
        cluster_group = "data/others/TCGA_IDHall_clusters.csv",
        kegg_pathway_genes_directory = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes_directory,
        genes_highExpressed = "data/others/transformed_expression_data/TCGA_genes_highExpressed.csv",

        #normalized_expression_data_expressionFiltered = rules.remove_low_expressed_genes.output.normalized_expression_data_expressionFiltered,

    output:
        DEgene_figures = directory("data/figures/TCGA_DEgene_figures_cluster_vs_others")
        
    conda:
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/extract_DEgenes_cluster_vs_others.R"
