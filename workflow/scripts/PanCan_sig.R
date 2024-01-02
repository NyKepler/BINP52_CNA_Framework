#!/usr/bin/env Rscript

# Title: PanCan_sig.R
# Author: Guyuan TANG
# Date: 2024/1/2

# Description: this script is used for the Pan-Cancer signature validation on each sample group.

library(tidyverse)
library(CINSignatureQuantification)

# specify the output directory
outdir <- snakemake@params[['outdir']]

# extract the required values from the snakemake workflow
groupID <- snakemake@params[['groupID']] # the group name
binsize <- snakemake@params[['binsize']] # the binsize for the selected group
sample_names <- snakemake@params[['sample_names']] # the list containing the sample names
in_dir <- snakemake@params[['indir']] # the directory to store the segment files

# generate the required sample table
sample_df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(sample_df) <- c('chromosome', 'start', 'end', 'segVal', 'sample')
## fill in the sample information
for (i in sample_names) {
  tab_df <- read.table(paste0(in_dir, i, '/05_absolute_CN/', i, '_', binsize, 'kb_seg.tsv'), sep = '\t', header = 1) %>% select('chromosome','start','end','segVal')
  tab_df['sample'] <- i
  ### append the table to the final sample table
  sample_df <- rbind(sample_df, tab_df)
}

# run the validation
PanCan_sig <- quantifyCNSignatures(sample_df)

# output the matrix files, matrix object(RDS) and a simple heatmap
output_prefix <- paste0(outdir, groupID, '/PanSig/', groupID, '_PanSig.', binsize, 'kb')
# matrix file
write.table(PanCan_sig@featFitting[["sampleByComponent"]], file = paste0(output_prefix, '.matrix.txt'))
# matrix object
saveRDS(PanCan_sig@featFitting[["sampleByComponent"]], file = paste0(output_prefix, '.matrix.rds'), compress = FALSE)
# object with full information (including activities, weights)
saveRDS(PanCan_sig, file = paste0(output_prefix, '.signature.rds'), compress = FALSE)
# simple heatmaps
pdf(paste0(output_prefix,'.activity.pdf'))
plotActivities(PanCan_sig)
dev.off()
pdf(paste0(output_prefix,'.component.pdf'))
plotSampleByComponent(PanCan_sig)
dev.off()
