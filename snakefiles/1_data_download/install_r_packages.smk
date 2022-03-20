
#######################################################

"""
            install R packages which are not available in the bionconda channel
"""

#######################################################

rule install_r_packages:
    output:
        r_package_status = ("data/others/r_package_status.csv")
    
    conda: 
        "../../envs/r.yaml"

    script:
        "../../scripts/r_modules/install_r_packages.R"
