
"""
    GO and KEGG enrichment analysis
"""

rule gene_enrichment_analysis:
    input:
        DEgene_results_dir = rules.performe_differential_gene_expression_analysis_cluster_vs_others.output.results_TCGA_IDHall_cluster_vs_others_dir,
        migration_genes = rules.extract_migration_markers.output.migration_genes
    output: 
        test = "data/test.png"
    conda:
        "../../envs/r.yaml"
    script:
        "../../scripts/r_modules/gene_enrichment_analysis.R" 