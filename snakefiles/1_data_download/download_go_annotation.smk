
"""
    Rule for downloading the annotations of the Gene Ontology
"""

rule download_go_annotation:
    output:
        go_annotation = "data/gene_ontology/homosapiens_goa.gaf"
    
    shell:
        """
        wget -O {output.go_annotation}.gz http://geneontology.org/gene-associations/goa_human.gaf.gz &&
        gzip -d {output.go_annotation}.gz
        """