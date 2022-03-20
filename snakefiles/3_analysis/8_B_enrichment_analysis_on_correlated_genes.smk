
"""
    Enrichment analysis on the correlated genes
"""

rule enrichment_analysis_on_correlated_genes:
    input:
        correlation_results = rules.calculate_gene_correlations.output.correlation_results
    output:
        dir_correlated_gene_enrichment_results = directory('data/{db}_{data_type}_correlated_gene_enrichment_results')
    params:
        utils_R = "scripts/r_modules/utils.R",
        ensembl_id = "ensembl_id.csv"  
    conda:
        "../../envs/r.yaml"
    script:
        "../../scripts/r_modules/enrichment_analysis_on_correlated_genes.R" 

