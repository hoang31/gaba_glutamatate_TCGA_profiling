
"""
    Analysis of the migration markers
"""

rule extract_migration_markers:
    input:
        migration_go_id = "migration_go_id.csv",
        go_annotation = rules.download_go_annotation.output.go_annotation

    output:
        migration_genes = "data/migration/migration_markers.csv"

    conda:
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/extract_migration_markers.R"