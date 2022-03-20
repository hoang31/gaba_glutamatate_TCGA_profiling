
######################################################
"""
            This rule is used for cleaning the data extracted from the CGGA database
"""
#######################################################

rule CGGA_cleaning:
    input:
        expression_data_CGGA1 = rules.download_CGGA1_data.output.CGGA1_expression,
        expression_data_CGGA2 = rules.download_CGGA2_data.output.CGGA2_expression,

        clinical_data_CGGA1 = rules.download_CGGA1_data.output.CGGA1_clinical,
        clinical_data_CGGA2 = rules.download_CGGA2_data.output.CGGA2_clinical,

    output:
        CGGA_expression_data = "data/expression_data/CGGA_expression_data.csv",
        CGGA_clinical_data = "data/clinical_data/CGGA_clinical_data.csv",

    script:
        "../../scripts/python_modules/CGGA_cleaning.py"
