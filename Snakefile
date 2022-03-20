
########################################
## Author Information
########################################


## author
__author__ = "Hoang Dong Nguyen"
__email__ = "hoang.dong.nguyen@usherbrooke.ca"


########################################
## Name of the project
########################################

"""
    Glioma characterization of GABA and Glutamate pathways
"""

########################################
## Import libraries and files
########################################


configfile: "config.json"
import os


########################################
## Paramaters
########################################


## number of threads which will be used
nb_threads = 2


########################################
## Wildcards
########################################


## wildcard for the TCGA, CGGA and other databases 
databases = [
    "TCGA"
]

data_type = [
    "IDHall",
    "IDHwt",
    "IDHmut"
]

data_type = [
    "IDHall"
]

########################################
## Import snakemake rules
########################################


#####
## Snakefiles for downloading
#####


include: 'snakefiles/1_data_download/download_CGGA1_data.smk'
include: 'snakefiles/1_data_download/download_CGGA2_data.smk'
include: 'snakefiles/1_data_download/download_IGAP_data.smk'
include: 'snakefiles/1_data_download/install_r_packages.smk'
include: 'snakefiles/1_data_download/download_go_annotation.smk'
include: 'snakefiles/1_data_download/download_immune_data_TCGA.smk'
include: 'snakefiles/1_data_download/download_immune_response_genes.smk'


#####
##  Snakefiles for cleaning the data
#####


include: 'snakefiles/2_data_cleaning/CGGA_cleaning.smk'
include: 'snakefiles/2_data_cleaning/IGAP_cleaning.smk'
include: 'snakefiles/2_data_cleaning/TCGA_cleaning.smk'


#####
## Snakefiles for analysis
#####


include: 'snakefiles/3_analysis/1_A_search_KEGG_PATHWAY_genes.smk'
include: 'snakefiles/3_analysis/1_B_estimate_purity.smk'
include: 'snakefiles/3_analysis/1_C_compare_KEGG_genes.smk'
include: 'snakefiles/3_analysis/2_normalize_expression_count.smk'
include: 'snakefiles/3_analysis/3_remove_low_expressed_genes.smk'
#include: 'snakefiles/3_analysis/4_remove_homogeneous_expressed_genes.smk'
include: 'snakefiles/3_analysis/5_generate_heatmap.smk'
include: 'snakefiles/3_analysis/6_A_performe_differential_gene_expression_analysis.smk'
include: 'snakefiles/3_analysis/6_B_calculate_survival.smk'
include: 'snakefiles/3_analysis/6_C_find_most_expressed_genes.smk'
include: 'snakefiles/3_analysis/6_D_analyse_epidemiological_data.smk'
include: 'snakefiles/3_analysis/6_E_analyse_epidemiological_data_subgroup.smk'
include: 'snakefiles/3_analysis/6_F_migration_markers.smk'
include: 'snakefiles/3_analysis/6_G_alteration_analysis.smk'
include: 'snakefiles/3_analysis/6_H_gene_enrichment_analysis.smk'
include: 'snakefiles/3_analysis/6_I_alteration_visualization.smk'
include: 'snakefiles/3_analysis/6_J_generate_expression_profile.smk'
include: 'snakefiles/3_analysis/6_K_analyze_methylation_data.smk'
include: 'snakefiles/3_analysis/6_L_generate_methylation_figures.smk'
include: 'snakefiles/3_analysis/6_M_generate_correlation_methylation_expression_figures.smk'
include: 'snakefiles/3_analysis/7_A_extract_DEgenes.smk'
include: 'snakefiles/3_analysis/7_B_extract_DEgenes_cluster_vs_others.smk'
include: 'snakefiles/3_analysis/8_A_calculate_gene_correlations.smk'
include: 'snakefiles/3_analysis/8_B_enrichment_analysis_on_correlated_genes.smk'
include: 'snakefiles/3_analysis/8_C_create_GO_subgroups.smk'
include: 'snakefiles/3_analysis/8_D_generate_figures_goterm_group.smk'
include: 'snakefiles/3_analysis/9_A_infer_immune_type_composition.smk'
include: 'snakefiles/3_analysis/9_B_pca_on_immune_data.smk'
include: 'snakefiles/3_analysis/10_A_validation_CGGA.smk'
include: 'snakefiles/3_analysis/generate_figures.smk'
include: 'snakefiles/3_analysis/generate_tables.smk'


#####
## snakefiles for the creation of the interaction network
#####


include: 'snakefiles/4_interaction_network/1_download_interaction_data.smk'
include: 'snakefiles/4_interaction_network/2_A_get_edge_information.smk'
include: 'snakefiles/4_interaction_network/2_B_get_node_information.smk'
include: 'snakefiles/4_interaction_network/3_A_IDH_gene_correlation.smk'
include: 'snakefiles/4_interaction_network/3_B_deseq2_miRNA.smk'
include: 'snakefiles/4_interaction_network/3_C_get_miRNA_interactions.smk'


########################################
## Rule all
########################################

rule all:
    input:
        ## generate the heatmap
        heatmap_IDHall = expand(
            os.path.join(
                "data/figures/",
                "{db}_heatmap_IDHall.svg"
            ),
            db = databases,
        ),

        ## extract_DEgenes rule
        upset_plot_DEgenes = expand(
            "data/figures/{db}_{data_type}_upset_plot_DEgenes.svg",
            db = databases,
            data_type = data_type
        ),

        volcano_plot_directory = expand(
            "data/figures/{db}_{data_type}_volcano_plot",
            db = databases,
            data_type = data_type
        ),

        expression_profile_plot_directory = expand(
                "data/figures/{db}_{data_type}_expression_profiles",
            db = databases,
            data_type = data_type    
        ),
        
        ## calculate_survival rule
        survival_curve_dir = expand(
            "data/figures/{db}_{data_type}_survival_curves",
            db = databases,
            data_type = data_type
        ),

        ## find_most_expressed_genes
        upset_plot = expand(
            "data/figures/{db}_{data_type}_upsetplot_topExpessedGenes.svg",
            db = databases,
            data_type = data_type
        ),

        ## generated figures
        generated_figures = expand(
            "data/figures/{db}_{data_type}_alteration_figures/",
            db = databases,
            data_type = data_type
        ),

        epidemio_dir = expand(
            "data/{db}_{data_type}_epidemiology",
            db = databases,
            data_type = data_type            
        ),

        ## analyse_epidemiologival_data_subgroup
        subgroup_epidemiological_data_all = expand(
            "data/{db}_{data_type}_epidemiology_subgroup/all_data.csv",
            db = databases,
            data_type = data_type
        ),

        ## rule 6_F_migration_markers
        migration_genes = "data/migration/migration_markers.csv",

        ## rule generate_figures_goterm_group
        #enrichment_figure_dir = directory('data/figures/TCGA_IDHall_enrichment'),

        ## rule 6_A_performe_differential_gene_analysis
        #results_TCGA_IDHall_cluster_vs_others_dir = ("data/DEGene_analysis_cluster_vs_others"),

        ## rule alteration_analysis
        alteration_visualization = expand(
            "data/figures/{db}_alteration_visualization.svg", 
            db = databases
        ),

        ### rule calculate_gene_correlations
        #correlation_results = expand(
        #    "data/others/{db}_{data_type}_correlation_results.csv",
        #    db = databases,
        #    data_type = data_type
        #),

        ### rule enrichment_analysis_on_correlated_genes
        #dir_correlated_gene_enrichment_results = expand(
        #    'data/{db}_{data_type}_correlated_gene_enrichment_results',
        #    db = databases,
        #    data_type = data_type
        #),

        ## rule create_GO_subgroups
        #go_groups = 'data/others/go_groups.csv',

        ### rule generate_figures_goterm_groups
        #enrichment_figure_dir = 'data/figures/TCGA_IDHall_enrichment',

        ## rule infer_immune_type_composition
        immune_composition_figures_dir = 'data/figures/TCGA_IDHall_immune_composition_figures',

        ##rule pca_on_immune_data
        pca_figure_dir = "data/figures/TCGA_IDHall_immune_pca",

        ## rule get_information_node
        information_node = config["path"]["information_node"],

        ## rule IDH_gene_correlation,
        #IDH_gene_correlation_dir = "data/IDH_gene_correlation",

        ## rule compare_KEGG_genes 
        upset_plot_genes_of_interest = "data/figures/upset_plot_genes_of_interest.svg",

        ## rule generate_tables
        descriptive_table_dir = "data/figures/TCGA_IDHall_descriptive_tables",

        ## rule generate_expression_profile
        neurotransmission_expression_profile_dir = "data/figures/TCGA_IDHall_neurotransmission_expression_profile",

        ## rule generate_methylation_figures
        methylation_figure_dir = "data/figures/TCGA_IDHall_methylation",
        methylation_expression_data = "data/others/gene_of_interest_data/1_methylation_expression_data.csv",

        ## rule generate_correlation_methylation_expression_figure
        methylation_figure_correlation_dir = "data/figures/TCGA_IDHall_methylation_correlation",

        ## rule extract_DEgenes_cluster_vs_others
        DEgene_figures = "data/figures/TCGA_DEgene_figures_cluster_vs_others",

        ## rules deseq2_miRNA
        results_deseq2_miRNA = "data/DEGene_analysis_miRNA",

        ## rules get_miRNA_interactions
        miRNA_figure_dir = "data/figures/TCGA_IDHall_miRNA",

        ## download and clean the CGGA data
        #cgga_validation_dir = "data/figures/CGGA_validation"



        #generate_figures2 = "data/generate_figures2"
        #test = "data/test.png"
