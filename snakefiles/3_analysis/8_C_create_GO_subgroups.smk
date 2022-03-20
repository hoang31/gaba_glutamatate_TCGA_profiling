
"""
    Create subgroup of Gene ontology terns for making the biological interpretation easier
"""

rule create_GO_subgroups:
    input:
        dir_correlated_gene_enrichment_results = ('data/TCGA_IDHall_correlated_gene_enrichment_results')
    output:
        go_groups = 'data/others/go_groups.csv'
    params:
        utils_R = "scripts/r_modules/utils.R"
    conda:
        "../../envs/r.yaml"
    script:
        "../../scripts/r_modules/create_GO_subgroups.R" 