
#######################################################
"""
            Download data from Ivy Glioblastoma Atlas Project (IGAP)
"""
#######################################################

rule download_IGAP_data:
    output:
        IGAP_expression = os.path.join(
            config["path"]["IGAP"], "IGAP_expression.txt"
        ),
        IGAP_clinical =  os.path.join(
            config["path"]["IGAP"], "IGAP_clinical.xls"
        )

    params:
        link_expression = config["download"]["IGAP_expression"],
        link_clinical = config["download"]["IGAP_clinical"],


    shell:
        """
        wget -O {output.IGAP_clinical} \"{params.link_clinical}\"
        """
