#!/usr/bin/env Rscript

# Title: PanCan_sig.R
# Author: Guyuan TANG
# Date: 2024/1/2

# Description: this script is used for the Pan-Cancer signature validation on each sample group.

library(tidyverse)
library(CINSignatureQuantification)
library(lsa)
library(openxlsx)

# specify the output directory
outdir <- snakemake@params[['outdir']]

# extract the required values from the snakemake workflow
sample_info <- snakemake@params[['sample_info']] # the sample information tsv file
sample_info <- read.table(sample_info, sep = '\t', header = 1)
sample_names <- sample_info$Sample
in_dir <- snakemake@params[['indir']] # the directory to store the segment files

# generate the required sample table
sample_df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(sample_df) <- c('chromosome', 'start', 'end', 'segVal', 'sample')
## fill in the sample information
for (i in sample_names) {
  binsize <- sample_info[which(sample_info$Sample==i), 'Binsize']
  tab_df <- read.table(paste0(in_dir, i, '/05_absolute_CN/', i, '_', binsize, 'kb_seg.tsv'), sep = '\t', header = 1) %>% select('chromosome','start','end','segVal')
  tab_df['sample'] <- i
  ### append the table to the final sample table
  sample_df <- rbind(sample_df, tab_df)
}

# run the validation
PanCan_sig <- quantifyCNSignatures(sample_df)

# output the matrix files, matrix object(RDS) and a simple heatmap
output_prefix <- paste0(outdir,'PanCan_sig')
# matrix file
write.table(PanCan_sig@featFitting[["sampleByComponent"]], file = paste0(output_prefix, '.SCmatrix.txt'))
# matrix object
saveRDS(PanCan_sig@featFitting[["sampleByComponent"]], file = paste0(output_prefix, '.SCmatrix.rds'), compress = FALSE)
# object with full information (including activities, weights)
saveRDS(PanCan_sig, file = paste0(output_prefix, '.signature.rds'), compress = FALSE)
# simple heatmaps
pdf(paste0(output_prefix,'.activity.pdf'))
plotActivities(PanCan_sig)
dev.off()
pdf(paste0(output_prefix,'.component.pdf'))
plotSampleByComponent(PanCan_sig)
dev.off()

# apply cosine similarity to find the enriched signature for each sample
## load the signature definition table for PanCancer signatures
pan_mat_input <- snakemake@params[['def_SC']]
pan_cancer_mat <- read.xlsx(pan_mat_input, colNames = TRUE, rowNames=TRUE)
pan_cancer_mat <- t(pan_cancer_mat)
pancan_sig <- colnames(pan_cancer_mat)

# calculation
SCmatrix <- PanCan_sig@featFitting[["sampleByComponent"]]
SCmatrix <- t(SCmatrix)
samples <- colnames(SCmatrix)
cos_mat <- cbind(pan_cancer_mat, SCmatrix)
cos_mat <- cosine(cos_mat)

# extract the required sample-by-signature similarity matrix
SSmatrix <- cos_mat[samples,]
SSmatrix <- as.data.frame(SSmatrix) %>% select(all_of(pancan_sig))

# find the most similar signature for each sample
sample_by_signature <- apply(SSmatrix, 1, function(t) colnames(SSmatrix)[which.max(t)])
sample_by_signature <- as.list(sample_by_signature)
SSmatrix['enrich'] <- NA
for (i in samples) {
  sigID <- sample_by_signature[[i]]
  SSmatrix[i, 'enrich'] <- sigID
}

# save the similarity SSmatrix
SSmatrix['sample'] <- rownames(SSmatrix)
SSmatrix <- SSmatrix %>% select('sample', 'enrich', 'CX1':'CX17')
row.names(SSmatrix) <- 1:nrow(SSmatrix)
write.table(SSmatrix, file = paste0(output_prefix,'.SSmatrix.tsv'), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
saveRDS(SSmatrix, file = paste0(output_prefix, '.SSmatrix.rds'), compress = FALSE)
