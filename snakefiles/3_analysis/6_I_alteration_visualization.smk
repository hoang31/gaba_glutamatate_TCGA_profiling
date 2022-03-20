

"""
    Add some known alerations for each cluster in the clustering
"""


rule alteration_visualization:
    input:
        kegg_pathway_genes = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes,        
        gene_id = rules.remove_low_expressed_genes.output.genes_highExpressed,
        expression_data = rules.normalize_expression_count.output.TCGA_expression_normalized,
        TCGA_clinical_cancer = "data/clinical_data/TCGA_clinical_data_cancer",
        TCGA_clinical_normal = "data/clinical_data/TCGA_clinical_data_normal",
        TCGA_clinical_EXTERNAL = "data/clinical_data/TCGA_clinical_EXTERNAL.csv",
        epidemiological_data = "data/TCGA_IDHall_epidemiology_subgroup/all_data.csv"
    output:
        alteration_visualization = "data/figures/{db}_alteration_visualization.svg",
        alteration_visualization_legend = "data/figures/{db}_alteration_visualization_legend.svg"
    params:
        heatmap_function = "scripts/r_modules/heatmap_generation.R",
        utils_R = "scripts/r_modules/utils.R",
        clustering_metrics = "scripts/r_modules/clustering_metrics.R",
    conda:
        "../../envs/r.yaml"
    script:
        "../../scripts/r_modules/alteration_visualization.R"    
