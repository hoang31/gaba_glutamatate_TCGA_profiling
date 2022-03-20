
######################################################
"""
            This rule is used for cleaning the data extracted from the CGGA database
"""
#######################################################

rule IGAP_cleaning:
    input:
        expression_data_IGAP = rules.download_IGAP_data.output.IGAP_expression,
        clinical_data_IGAP = rules.download_IGAP_data.output.IGAP_clinical,


    output:
        IGAP_expression_data = "data/expression_data/IGAP_expression_data.csv",
        IGAP_clinical_data = "data/clinical_data/IGAP_clinical_data.csv",

    script:
        "../../scripts/python_modules/IGAP_cleaning.py"
