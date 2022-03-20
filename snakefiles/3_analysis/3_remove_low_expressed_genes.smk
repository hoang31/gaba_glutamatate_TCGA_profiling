
######################################################
"""
            It permits to filter the expression data extracting GABA and GLU genes and removing the genes which are lowly expressed
"""
#######################################################

rule remove_low_expressed_genes:
    input:
        normalized_expression_data = rules.normalize_expression_count.output.TCGA_expression_normalized,
        kegg_pathway_genes = rules.search_KEGG_PATHWAY_genes.output.kegg_pathway_genes,

    output:
        genes_highExpressed = "data/others/transformed_expression_data/{db}_genes_highExpressed.csv",
        normalized_expression_data_expressionFiltered = os.path.join(
            config["path"]["expression_data"],
            "{db}_normalized_expression_data_expressionFiltered.csv"
        ),
        normalized_expression_data_sample_data = os.path.join(
            config["path"]["clinical_data"],
            "{db}_normalized_expression_data_sample_data.csv"
        ),
        
    params:
        gene_id = "ensembl_id.csv",
        expression_filter = 10,
        utils_R = "scripts/r_modules/utils.R",

    conda:
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/remove_low_expressed_genes.R"
