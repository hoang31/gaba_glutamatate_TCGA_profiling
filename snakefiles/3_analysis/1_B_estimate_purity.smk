#######################################################

"""
        Estimate the purity for each sample from the FPKM expression data
"""

#######################################################

rule estimate_purity:
    input:
        r_package_status = rules.install_r_packages.output.r_package_status,
        fpkm_expression_data = "data/expression_data/{db}_expression_data_fpkm_cancer",
        normal_fpkm_expression_data = "data/expression_data/TCGA_expression_data_fpkm_normal"

    output:
        purity_score_tumor = "data/others/{db}_purity_score_tumor.csv",
        purity_score_normal = "data/others/{db}_purity_score_normal.csv",

    params:
        utils_R = "scripts/r_modules/utils.R",
        ensembl_id = "ensembl_id.csv"

    conda: 
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/estimate_purity.R"
