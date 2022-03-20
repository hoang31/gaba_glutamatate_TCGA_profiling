
"""
    Download the interaction data from databases
"""

rule download_stringdb_protein_protein:
    output:
        stringdb_protein_protein = config["path"]["stringdb_protein_protein"]

    params:
        link = config["download"]["stringdb_protein_protein"]

    shell:
        "wget -O {output.stringdb_protein_protein}.gz {params.link}"
        " && gunzip {output.stringdb_protein_protein}.gz"


rule download_stringdb_protein_information:
    output:
        stringdb_protein_information = config["path"]["stringdb_protein_information"]

    params:
        link = config["download"]["stringdb_protein_information"]

    shell:
        "wget -O {output.stringdb_protein_information}.gz {params.link}"
        " && gunzip {output.stringdb_protein_information}.gz"


rule download_rnainter_rna_rna:
    output:
        rnainter_rna_rna = config["path"]["rnainter_rna_rna"]

    params:
        link = config["download"]["rnainter_rna_rna"]

    shell:
        "wget -O {output.rnainter_rna_rna}.zip {params.link}"
        " && unzip {output.rnainter_rna_rna}.zip"
        " && mv RNA-RNA.txt {output.rnainter_rna_rna}"


rule download_rnainter_rna_protein:
    output:
        rnainter_rna_protein = config["path"]["rnainter_rna_protein"]

    params:
        link = config["download"]["rnainter_rna_protein"]

    shell:
        "wget -O {output.rnainter_rna_protein}.zip {params.link}"
        " && unzip {output.rnainter_rna_protein}.zip"
        " && mv RNA-Protein.txt {output.rnainter_rna_protein}"

rule filter_rna_protein_interaction_data:
    input:
        rnainter_rna_protein = rules.download_rnainter_rna_protein.output.rnainter_rna_protein
    output:
        rnainter_rna_protein_filtered = config["path"]["rnainter_rna_protein_filtered"]
    params:
        confidence_score_cutoff = 0.5
    conda:
        "../../envs/r.yaml"
    script:
        "../../scripts/r_modules/filter_rna_protein_interaction_data.R"


rule download_rnainter_rna_compound:
    output:
        rnainter_rna_compound = config["path"]["rnainter_rna_compound"]

    params:
        link = config["download"]["rnainter_rna_compound"]

    shell:
        "wget -O {output.rnainter_rna_compound}.zip {params.link}"
        " && unzip {output.rnainter_rna_compound}.zip"
        " && mv RNAInter_interaction_RC.txt {output.rnainter_rna_compound}"