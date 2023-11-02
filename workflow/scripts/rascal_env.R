# this script is used for setting up the conda R environment with rascal package

library(devtools)

install_github("crukci-bioinformatics/rascal")

sink(file = "other_info/rascal_settle_info.txt")
print("build-up rascal package and environment")
sink(file = NULL)