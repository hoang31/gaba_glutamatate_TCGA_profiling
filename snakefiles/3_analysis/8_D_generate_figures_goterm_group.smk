

"""
    Generate figures for go term group visualization
"""

rule generate_figures_goterm_group:
    input:
        go_groups = 'data/others/go_groups.csv'

    output:
        enrichment_figure_dir = directory('data/figures/TCGA_IDHall_enrichment')
    params:
        utils_R = "scripts/r_modules/utils.R"
    conda:
        "../../envs/r.yaml"
    script:
        "../../scripts/r_modules/generate_figures_goterm_group.R" 