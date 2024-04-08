# Title: panConusig_env.R
# Author: Guyuan TANG
# Date: 2024-03-05

# Description: this is used to set up the environment with required packages to do the panConusig analysis.

# Battenburg and ASCAT.sc
BiocManager::install(c( "minfi", "conumee", "Rsamtools", "igordot/copynumber"), update = TRUE, ask = FALSE) # "copynumber"

devtools::install_github('VanLoo-lab/ascat/ASCAT', upgrade = TRUE)
devtools::install_github("Wedge-Oxford/battenberg", upgrade = TRUE)
devtools::install_github("VanLoo-lab/ASCAT.sc", build_opts = c("--no-build-vignettes"), upgrade = TRUE)
devtools::install_github("UCL-Research-Department-of-Pathology/panConusig")

# test and print
library(Battenberg)
library(ASCAT.sc)
library(panConusig)
