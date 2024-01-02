#!/usr/bin/env Rscript
# this script is used for setting up the conda R environment required for running the pan-cancer signature validation.

library(devtools)

sink(file = snakemake@output[[1]])
install_github("markowetzlab/CINSignatureQuantification", build_vignettes = TRUE, dependencies = TRUE)
library(CINSignatureQuantification)
print("build-up required environment for Pan-Cancer signatures validation!")
sink(file = NULL)