
"""
    Download the genes related to the Gene ontoly term of immune response (go id GO:0002376)
"""

rule download_immune_response_genes:
    output:
        immune_response_genes = "data/others/immune_response_genes.txt"
    
    params:
        link = config["download"]['immune_response_genes']
    
    shell:
        "wget -O {output.immune_response_genes} '{params.link}'"
