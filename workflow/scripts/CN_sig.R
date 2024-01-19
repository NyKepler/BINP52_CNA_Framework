#!/usr/bin/env Rscript

# Title: CN_sig.R
# Author: Guyuan TANG
# Date: 2023/12/29

# Description: this script is used for signature validation on each sample group.

library(NMF)
library(flexmix)
library(YAPSA)
library(tidyverse)
library(lsa)

main_script = snakemake@input[['main']]
helper_script = snakemake@input[['helper']]

# load the functions
source(main_script)
source(helper_script)

# specify the output directory
outdir <- snakemake@params[['outdir']]

# extract the required values from the snakemake workflow
sample_info <- snakemake@params[['sample_info']] # the sample information tsv file
sample_info <- read.table(sample_info, sep = '\t', header = 1)
sample_names <- sample_info$Sample
in_dir <- snakemake@params[['indir']] # the directory to store the segment files

# generate the required input dataframes list
segment_list <- list()
for (i in sample_names) {
  binsize <- sample_info[which(sample_info$Sample==i),'Binsize']
  tab_df <- read.table(paste0(in_dir, i, '/05_absolute_CN/', i, '_', binsize, 'kb_seg.tsv'), sep = '\t', header = 1) %>% select('chromosome','start','end','segVal')
  a <- length(segment_list) + 1
  segment_list[[a]] <- tab_df
}
names(segment_list) <- sample_names

# validate the signatures in the samples
CN_features <- extractCopynumberFeatures(segment_list)
sample_by_component <- generateSampleByComponentMatrix(CN_features)

# output the matrix files, matrix object(RDS) and a simple heatmap
output_prefix <- paste0(outdir,'CN_sig')
# the same output for the sample-by-component matrix should also be saved for future analysis
write.table(sample_by_component, file = paste0(output_prefix, '.SCmatrix.txt'))
saveRDS(sample_by_component, file = paste0(output_prefix, '.SCmatrix.rds'), compress = FALSE)
pdf(paste0(output_prefix, 'heatmapSC.pdf'))
heatmap(sample_by_component)
dev.off()


# apply cosine similarity to find the enriched signature for each sample
## load the signature definition table for Brenton's CN signatures
feat_sig_mat <- snakemake@params[['def_SC']]
CN_sig <- paste0('s',c(1:7))

## load the Sample-by-Component matrix generated above
SCmatrix <- t(sample_by_component)
samples <- colnames(SCmatrix)

## calculate the cosine similarity
cos_mat <- cbind(feat_sig_mat, SCmatrix)
cos_mat <- cosine(cos_mat)

## extract the required signature similarity matrix
SSmatrix <- cos_mat[samples,]
SSmatrix <- as.data.frame(SSmatrix) %>% select(all_of(CN_sig))

## find the most similar signature for each sample
sample_by_signature <- apply(SSmatrix, 1, function(t) colnames(SSmatrix)[which.max(t)]) 
sample_by_signature <- as.list(sample_by_signature)
SSmatrix['enrich'] <- NA
for (i in samples) {
  sigID <- sample_by_signature[[i]]
  SSmatrix[i, 'enrich'] <- sigID
}

## save the SSmatrix
SSmatrix['sample'] <- rownames(SSmatrix)
SSmatrix <- SSmatrix %>% select('sample', 'enrich', 's1':'s7')
row.names(SSmatrix) <- 1:nrow(SSmatrix)
write.table(SSmatrix, file = paste0(output_prefix,'.SSmatrix.tsv'), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
saveRDS(SSmatrix, file = paste0(output_prefix,'.SSmatrix.rds'), compress = FALSE)



