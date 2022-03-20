
######################################################
"""
            COPY NUMBER AND SINGLE NUCLEOTIDE VARIATION analysis
"""
#######################################################

rule alteration_analysis:
    input:
        cnv_gbm = "data/copy_number_variation_data/GBM.focal_score_by_genes.txt",
        cnv_lgg = "data/copy_number_variation_data/LGG.focal_score_by_genes.txt",
        cnv_id = "data/copy_number_variation_data/gdc_sample_sheet.2020-05-06.tsv",
        # snv_data = "data/single_nucleotide_variation_data/",
        gene_id = rules.remove_low_expressed_genes.output.genes_highExpressed,
        # expression_data = "data/{db}_expression_normalized.csv",
        clinical_data = "data/clinical_data/TCGA_clinical_data_cancer",
        # DEgene_results_dir = rules.performe_differential_gene_expression_analysis.output.DEgene_results_dir,

    output:
        alteration_data = "data/others/{db}_{data_type}_alteration_data.csv"

    params:
        utils_R = "scripts/r_modules/utils.R",

    conda:
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/alteration_analysis.R"
