#!/usr/bin/env Rscript

# Title: CN_sig.R
# Author: Guyuan TANG
# Date: 2023/12/29

# Description: this script is used for signature validation on each sample group.

library(NMF)
library(flexmix)
library(YAPSA)
library(tidyverse)

main_script = snakemake@input[['main']]
helper_script = snakemake@input[['helper']]

# load the functions
source(main_script)
source(helper_script)

# specify the output directory
outdir <- snakemake@params[['outdir']]

# extract the required values from the snakemake workflow
groupID <- snakemake@params[['groupID']] # the group name
binsize <- snakemake@params[['binsize']] # the binsize for the selected group
sample_names <- snakemake@params[['sample_names']] # the list containing the sample names
in_dir <- snakemake@params[['indir']] # the directory to store the segment files

# generate the required input dataframes list
segment_list <- list()
for (i in sample_names) {
  tab_df <- read.table(paste0(in_dir, i, '/05_absolute_CN/', i, '_', binsize, 'kb_seg.tsv'), sep = '\t', header = 1) %>% select('chromosome','start','end','segVal')
  a <- length(segment_list) + 1
  segment_list[[a]] <- tab_df
}
names(segment_list) <- sample_names

# validate the signatures in the samples
CN_features <- extractCopynumberFeatures(segment_list)
sample_by_component <- generateSampleByComponentMatrix(CN_features)
output_matrix <- quantifySignatures(sample_by_component)

# output the matrix files, matrix object(RDS) and a simple heatmap
output_prefix <- paste0(outdir, groupID, '/CNsig/', groupID, '_CNsig.', binsize, 'kb')
# matrix file
write.table(output_matrix, file = paste0(output_prefix, '.matrix.txt'))
# matrix object
saveRDS(output_matrix, file = paste0(output_prefix, '.matrix.rds'), compress = FALSE)
# simple heatmap
pdf(paste0(output_prefix,'.heatmapSS.pdf'))
heatmap(output_matrix)
dev.off()
# the same output for the sample-by-component matrix should also be saved for future analysis
write.table(sample_by_component, file = paste0(output_prefix, '.SCmatrix.txt'))
saveRDS(sample_by_component, file = paste0(output_prefix, '.SCmatrix.rds'), compress = FALSE)
pdf(paste0(output_prefix, 'heatmapSC.pdf'))
heatmap(sample_by_component)
dev.off()


