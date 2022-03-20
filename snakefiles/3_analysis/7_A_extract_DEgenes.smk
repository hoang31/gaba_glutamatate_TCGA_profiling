
######################################################

"""
            Analyse the differential expressed genes
"""

#######################################################


rule extract_DEgenes:
    input:
        DEgene_results_dir = rules.performe_differential_gene_expression_analysis.output.DEgene_results_dir,
        cluster_group = "data/others/{db}_{data_type}_clusters.csv",
        normalized_expression_data_expressionFiltered = rules.remove_low_expressed_genes.output.normalized_expression_data_expressionFiltered,
        kegg_pathway_genes_directory = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes_directory,
        normal_fpkm_expression_data = "data/expression_data/TCGA_expression_data_fpkm_normal"

    output:
        upset_plot_DEgenes = (
            "data/figures/{db}_{data_type}_upset_plot_DEgenes.svg"
        ),

        volcano_plot_directory = directory(
            "data/figures/{db}_{data_type}_volcano_plot"
        ),

        expression_profile_plot_directory = directory(
            "data/figures/{db}_{data_type}_expression_profiles"
        ),

        DE_genes = "data/others/{db}_{data_type}_DE_genes.txt",
        DE_genes_in_common = "data/others/{db}_{data_type}_DE_genes_in_common.txt",

    params:
        wcards = "{data_type}",
        utils_R = "scripts/r_modules/utils.R",
        ensembl_id = "ensembl_id.csv"
        
    conda:
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/extract_DEgenes.R"
