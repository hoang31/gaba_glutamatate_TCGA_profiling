
#######################################################
"""
            Download data from Chinese Glioma Genome Atlas (CGGA)
"""
#######################################################

rule download_CGGA2_data:
    output:
        CGGA2_expression = os.path.join(
            config["path"]["CGGA2"], "CGGA2_expression.txt"
        ),
        CGGA2_clinical =  os.path.join(
            config["path"]["CGGA2"], "CGGA2_clinical.txt"
        )

    params:
        link_expression = config["download"]["CGGA2_expression"],
        link_clinical = config["download"]["CGGA2_clinical"],


    shell:
        """
        wget -O {output.CGGA2_expression}.zip \"{params.link_expression}\" &&
        zip -d {output.CGGA2_expression}.zip __MACOSX/\* &&
        unzip -c {output.CGGA2_expression}.zip | tail -n +3 > {output.CGGA2_expression} &&

        wget -O {output.CGGA2_clinical}.zip \"{params.link_clinical}\" &&
        zip -d {output.CGGA2_clinical}.zip __MACOSX/\* &&
        unzip -c {output.CGGA2_clinical}.zip | tail -n +3 > {output.CGGA2_clinical} &&


        rm {output.CGGA2_expression}.zip &&
        rm {output.CGGA2_clinical}.zip
        """
