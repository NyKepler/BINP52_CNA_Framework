# this script is used for setting up the conda R environment with rascal package

library(devtools)

sink(file = snakemake@output[[1]])
install_github("crukci-bioinformatics/rascal")
print("build-up rascal package and environment")
sink(file = NULL)