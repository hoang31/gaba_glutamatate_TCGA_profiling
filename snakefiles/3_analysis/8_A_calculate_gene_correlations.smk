
"""
    calculate the correlations of the DE genes with the other genes
"""

rule calculate_gene_correlations:
    input:
        gene_of_interest = rules.extract_DEgenes.output.DE_genes_in_common,
        TCGA_expression_normalized = rules.normalize_expression_count.output.TCGA_expression_normalized,
        TCGA_deseq_res_dir = "data/DEGene_analysis/results_TCGA_IDHall"
    output:
        correlation_results = "data/others/{db}_{data_type}_correlation_results.csv"
    params:
        ensembl_id = "ensembl_id.csv"
    conda:
        "../../envs/r.yaml" 
    script:
        "../../scripts/r_modules/calculate_gene_correlations.R" 
