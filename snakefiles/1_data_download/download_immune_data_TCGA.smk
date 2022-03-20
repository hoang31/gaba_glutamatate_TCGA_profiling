
"""
        Download the immune data of the TCGA samples from the TIMER2 webserver
"""

rule download_immune_data_TCGA:
    output:
        TCGA_immune_data = "data/others/TCGA_immune_data.tsv"
    
    params:
        link = config["download"]["TCGA_immune_composition"]

    shell:
        "wget -O {output.TCGA_immune_data}.gz {params.link} && "
        "gzip -d {output.TCGA_immune_data}.gz "