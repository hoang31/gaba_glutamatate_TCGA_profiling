
######################################################
"""
            It permits to calculate the local coefficient of variation
"""
#######################################################

rule remove_homogeneous_expressed_genes:
    input:
        expression_data = "data/{db}_expression_normalized.csv",
        sample_data = rules.remove_low_expressed_genes.output.normalized_expression_data_sample_data,
        genes_highExpressed = rules.remove_low_expressed_genes.output.genes_highExpressed,

    output:
        genes_highExpressed_heteregenousExpressed = os.path.join(
            config["path"]["others"],        
            "{db}_genes_highExpressed_heteregenousExpressed.csv"
    )

    params:
        window_size = 50,
        percentile_filter_value = 0.50

    threads:
        nb_threads

    conda:
        "../../envs/py_env.yaml"

    script:
        "../../scripts/python_modules/remove_homogeneous_expressed_genes.py"
