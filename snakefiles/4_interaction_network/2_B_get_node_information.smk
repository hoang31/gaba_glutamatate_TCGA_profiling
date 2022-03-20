
"""
        Clean the interaction data that were downloaded
"""

rule get_node_information:
    input:
        information_edge = rules.get_edge_information.output.information_edge,
        immune_response_genes = rules.download_immune_response_genes.output.immune_response_genes,
        genes_highExpressed = "data/others/transformed_expression_data/TCGA_genes_highExpressed.csv"

    output:
        information_node = config["path"]["information_node"]
    
    params:
        kegg_pathway_dir = "data/others/kegg_pathway_genes_directory"

    conda:
        "../../envs/py_env.yaml"

    script:
        "../../scripts/python_modules/get_node_information.py"