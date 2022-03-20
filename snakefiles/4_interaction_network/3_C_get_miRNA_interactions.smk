
"""
    Get the interactions between neurotransmission related genes and the miRNAs
"""

rule get_miRNA_interactions:
    input:
        methylation_expression_data = rules.generate_correlation_methylation_expression_figures.output.methylation_expression_data,
        #neurotransmission_gene_of_interest = "data/others/transformed_expression_data/TCGA_genes_highExpressed.csv",
        edge_information = rules.get_edge_information.output.information_edge,
        cluster_group = "data/others/TCGA_IDHall_clusters.csv",
        results_deseq2_miRNA = rules.deseq2_miRNA.output.results_deseq2_miRNA,
        miRNA_expression_data_healthy = "data/expression_data/TCGA_miRNA_expression_healthy_samples.csv",

        miRNA_expression_data = "data/others/TCGA_deseq2_mir_normalized_expression.csv",
        #miRNA_expression_data = "data/expression_data/TCGA_miRNA_expression.csv",


    output:
        miRNA_figure_dir = directory("data/figures/TCGA_IDHall_miRNA")

    params:
        utils = "scripts/r_modules/utils.R",
        pval_cutoff = 0.001,
        logfoldchange_cutoff = 1,
        expression_cutoff = 1

    conda:
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/get_miRNA_interactions.R"
