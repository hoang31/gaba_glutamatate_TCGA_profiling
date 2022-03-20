
#########################################################
#########################################################

"""
        SURVIVAL ANALYSIS of each clusters
"""


#########################################################
#########################################################


rule find_most_expressed_genes:
    input: 
        cluster_group = "data/others/{db}_{data_type}_clusters.csv",
        expression_data = rules.normalize_expression_count.output.TCGA_expression_normalized,
        gene_id = rules.remove_low_expressed_genes.output.genes_highExpressed,

    output: 
       upset_plot = "data/figures/{db}_{data_type}_upsetplot_topExpessedGenes.svg"

    params: 
        wcard = "{data_type}",
        utils_R = "scripts/r_modules/utils.R",
        top_Expressed_Gene_threshold = 50

    conda:
        "../../envs/r.yaml"

    script: 
        "../../scripts/r_modules/find_most_expressed_genes.R"