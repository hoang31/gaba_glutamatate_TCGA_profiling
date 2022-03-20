
"""
            This rule will normalize the expression data using the function rlog()/vst() available in the DESEQ2 R packages.
"""

rule normalize_expression_count:
    input:
        CGGA_expression = "data/expression_data/CGGA_expression_data.csv",
        CGGA_clinical = "data/clinical_data/CGGA_clinical_data.csv",

        TCGA_expression = "data/expression_data/TCGA_expression_data_count_cancer",
        TCGA_clinical = "data/clinical_data/TCGA_clinical_data_cancer",
        TCGA_expression_normal = "data/expression_data/TCGA_expression_data_count_normal",
        TCGA_clinical_normal = "data/clinical_data/TCGA_clinical_data_normal",

        IGAP_expression = "data/expression_data/IGAP_expression_data.csv",
        IGAP_clinical = "data/clinical_data/IGAP_clinical_data.csv",

        GTEX_expression = "data/expression_data/GTEX_expression_data.csv",
        GTEX_clinical = "data/clinical_data/GTEX_clinical_data.csv",

    output:
        TCGA_expression_normalized = "data/transformed_expression_data/TCGA_expression_normalized.csv"

    params:
        utils_R = "scripts/r_modules/utils.R",
        TCGA_clinical_EXTERNAL = "data/clinical_data/TCGA_clinical_EXTERNAL.csv",

    conda:
        "../../envs/r_deseq2.yaml"

    script:
        "../../scripts/r_modules/normalize_expression_count.R"
