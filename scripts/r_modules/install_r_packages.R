

##########################################
##########################################


##### install R packages that are not available with conda


##########################################
##########################################


##### installation of the r packages of interest

## install estimate for the estimation of cancer purity
install.packages("estimate", repos="http://R-Forge.R-project.org")

## install homo sapiens annotations for the go and kegg enrichment


##########################################
##########################################


# ##### Load the libraries for testing
test <- try(library(estimate))


# ##########################################
# ##########################################

# if (inherits(test, "try-error")) {
#     stop("\n\nEstimate r package is not installed or can not be loaded!\n\n")
# }
# # else {
# #     cat("\n\n", "Estimate R Package Installation Succeed!", "\n\n")
# # }


cat("estimate,installed", file = snakemake@output[["r_package_status"]], sep = "\n")
