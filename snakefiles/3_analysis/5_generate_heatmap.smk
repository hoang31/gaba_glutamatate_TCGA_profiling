
######################################################
"""
            It permits to create heatmap using the expression data and the GABA or the GLUTA genes
"""
#######################################################

rule generate_heatmap:
    input:
        kegg_pathway_genes = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes,        
        gene_id = rules.remove_low_expressed_genes.output.genes_highExpressed,
        expression_data = rules.normalize_expression_count.output.TCGA_expression_normalized,
        TCGA_clinical_cancer = "data/clinical_data/TCGA_clinical_data_cancer",
        TCGA_clinical_normal = "data/clinical_data/TCGA_clinical_data_normal",
        TCGA_clinical_EXTERNAL = "data/clinical_data/TCGA_clinical_EXTERNAL.csv",
        kegg_pathway_genes_directory = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes_directory,


    output:
        ## heatmaps
        heatmap_IDHall= os.path.join(
            "data/figures/",
            "{db}_heatmap_IDHall.svg"
        ),
        # heatmap_IDHwt = os.path.join(
        #     "data/figures/",
        #     "{db}_heatmap_IDHwt.svg"
        # ),
        # heatmap_IDHmut = os.path.join(
        #     "data/figures/",
        #     "{db}_heatmap_IDHmut.svg"
        # ),

        ## clusters groups
        cluster_IDHall = os.path.join(
            "data/others/",
            "{db}_IDHall_clusters.csv"
        ),

        #cluster_IDHall_subclusters = os.path.join(
        #    "data/others/",
        #    "{db}_IDHall_subclusters.csv"
        #),

        # cluster_IDHwt = os.path.join(
        #     "data/others/",
        #     "{db}_IDHwt_clusters.csv"
        # ),
        # cluster_IDHmut = os.path.join(
        #     "data/others/",
        #     "{db}_IDHmut_clusters.csv"
        # ),

        optimal_nb_clusters_IDHall = os.path.join(
            "data/figures/",
            "{db}_optimal_nb_clusters_IDHall.svg"
        ),

        # optimal_nb_clusters_IDHwt = os.path.join(
        #     "data/figures/",
        #     "{db}_optimal_nb_clusters_IDHwt.svg"
        # ),
        # optimal_nb_clusters_IDHmut = os.path.join(
        #     "data/figures/",
        #     "{db}_optimal_nb_clusters_IDHmut.svg"
        # ),
       
    params:
        heatmap_function = "scripts/r_modules/heatmap_generation.R",
        utils_R = "scripts/r_modules/utils.R",
        clustering_metrics = "scripts/r_modules/clustering_metrics.R",
        

    conda:
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/generate_heatmap.R"
