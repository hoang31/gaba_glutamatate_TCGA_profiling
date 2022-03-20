
"""
    Calculate the correlation between the IDH expression and other interesting genes (identified by network analysis)
"""

rule  IDH_gene_correlation:
    input:
        cluster_group = "data/others/TCGA_IDHall_clusters.csv",
        TCGA_expression_normalized = rules.normalize_expression_count.output.TCGA_expression_normalized,

    output:
        IDH_gene_correlation_dir = directory("data/IDH_gene_correlation")

    params:
        utils_R = "scripts/r_modules/utils.R"

    conda:
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/IDH_gene_correlation.R"