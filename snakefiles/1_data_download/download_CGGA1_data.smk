
#######################################################
"""
            Download data from Chinese Glioma Genome Atlas (CGGA)
"""
#######################################################

rule download_CGGA1_data:
    output:
        CGGA1_expression = os.path.join(
            config["path"]["CGGA1"], "CGGA1_expression.txt"
        ),
        CGGA1_clinical =  os.path.join(
            config["path"]["CGGA1"], "CGGA1_clinical.txt"
        )

    params:
        link_expression = config["download"]["CGGA1_expression"],
        link_clinical = config["download"]["CGGA1_clinical"],


    shell:
        """
        wget -O {output.CGGA1_expression}.zip \"{params.link_expression}\" &&
        zip -d {output.CGGA1_expression}.zip __MACOSX/\* &&
        unzip -c {output.CGGA1_expression}.zip | tail -n +3 > {output.CGGA1_expression} &&

        wget -O {output.CGGA1_clinical}.zip \"{params.link_clinical}\" &&
        zip -d {output.CGGA1_clinical}.zip __MACOSX/\* &&
        unzip -c {output.CGGA1_clinical}.zip | tail -n +3 > {output.CGGA1_clinical} &&


        rm {output.CGGA1_expression}.zip &&
        rm {output.CGGA1_clinical}.zip

        """
