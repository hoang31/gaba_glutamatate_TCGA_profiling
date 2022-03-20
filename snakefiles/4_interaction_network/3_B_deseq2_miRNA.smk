
"""
    Rule for the Differential miRNA gene expression using Deseq2
"""


rule deseq2_miRNA:
    input:
        edge_information = rules.get_edge_information.output.information_edge,
        miRNA_expression_data_rawcount = "data/expression_data/TCGA_miRNA_expression_rawcount.csv",
        cluster_group = "data/others/TCGA_IDHall_clusters.csv",

        mir_expression_healthy_dir = "data/expression_data/TCGA_mir_healthy_samples"

    output:
        results_deseq2_miRNA = directory("data/DEGene_analysis_miRNA"),
        deseq2_miRNA_normalized_expression = "data/others/TCGA_deseq2_mir_normalized_expression.csv"

    params:
        utils = "scripts/r_modules/utils.R",
        pval_cutoff = 0.001,
        logfoldchange_cutoff = 0

    conda:
        "../../envs/r_deseq2.yaml"

    script:
        "../../scripts/r_modules/deseq2_miRNA.R"
