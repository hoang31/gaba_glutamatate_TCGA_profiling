
"""
        Clean the interaction data that were downloaded
"""

rule get_edge_information:
    input:
        stringdb_protein_protein = rules.download_stringdb_protein_protein.output.stringdb_protein_protein,
        stringdb_protein_information = rules.download_stringdb_protein_information.output.stringdb_protein_information,

        rnainter_rna_rna = rules.download_rnainter_rna_rna.output.rnainter_rna_rna,
        rnainter_rna_protein_filtered = rules.filter_rna_protein_interaction_data.output.rnainter_rna_protein_filtered,
        rnainter_rna_compound = rules.download_rnainter_rna_compound.output.rnainter_rna_compound,

    output:
        information_edge = config["path"]["information_edge"]

    conda:
        "../../envs/py_env.yaml"

    script:
        "../../scripts/python_modules/get_edge_information.py"