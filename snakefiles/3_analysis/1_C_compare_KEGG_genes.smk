
######################################################
"""
            Compare all the genes extracted from the KEGG pathway databases
"""
#######################################################


rule compare_KEGG_genes:
    input:
        kegg_pathway_genes_directory = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes_directory,

    output:
        upset_plot_genes_of_interest = "data/figures/upset_plot_genes_of_interest.svg",

    
    conda:
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/compare_KEGG_genes.R"
