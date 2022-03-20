
######################################################

"""
            Differiential gene expression analysis using DESEQ2 R packages 
"""

#######################################################


rule performe_differential_gene_expression_analysis:
    input:
        TCGA_expression = "data/expression_data/TCGA_expression_data_count_cancer",
        TCGA_clinical = "data/clinical_data/TCGA_clinical_data_cancer",
        cluster_group = "data/others/TCGA_IDHall_clusters.csv",
        gene_of_interest = rules.remove_low_expressed_genes.output.genes_highExpressed,

    output:
        DEgene_results_dir = directory(
            "data/DEGene_analysis/results_{db}_{data_type}"
        )
       
    params:
        utils_R = "scripts/r_modules/utils.R",
        ensembl_id = "ensembl_id.csv"       

    conda:
        "../../envs/r_deseq2.yaml"

    script:
        "../../scripts/r_modules/performe_differential_gene_expression_analysis.R"



rule performe_differential_gene_expression_analysis_cluster_vs_others:
    input:
        TCGA_expression = "data/expression_data/TCGA_expression_data_count_cancer",
        TCGA_clinical = "data/clinical_data/TCGA_clinical_data_cancer",
        cluster_group = "data/others/TCGA_IDHall_clusters.csv",

    output:
        results_TCGA_IDHall_cluster_vs_others_dir = 
            directory("data/DEGene_analysis_cluster_vs_others")
       
    params:
        utils_R = "scripts/r_modules/utils.R",
        ensembl_id = "ensembl_id.csv"       

    conda:
        "../../envs/r_deseq2.yaml"

    script:
        "../../scripts/r_modules/performe_differential_gene_expression_analysis_cluster_vs_others.R"
