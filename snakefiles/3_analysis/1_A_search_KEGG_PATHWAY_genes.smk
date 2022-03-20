
######################################################
"""
            It permits to search the genes associated with the pathway using KEGG database
"""
#######################################################

rule search_KEGG_PATHWAY_genes:
    output:
        kegg_pathway_genes = "data/others/kegg_pathway_genes.txt",
        kegg_pathway_genes_directory = directory("data/others/kegg_pathway_genes_directory")
        
    script:
        "../../scripts/python_modules/search_genes_from_KEGGpathways.py"
