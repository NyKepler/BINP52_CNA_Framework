# Title: panConusig.R
# Author: Guyuan TANG
# Date: 2024/03/11

# Description: this is used to generate the sample-by-component matrix and the sample-by-signature matrix (cosine similarity) of panConusig, based on the allele-specific copy number profiles.

library(panConusig)
library(tidyverse)
library(lsa)

# load the copy number profiles for all the samples
input_file_list <- snakemake@input[['cna_profiles']]
in_mat <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(in_mat) <- c('sample', 'chr', 'startpos', 'endpos', 'nA', 'nB')
for (cna_profile in input_file_list) {
  df <- read.csv(cna_profile, sep = '\t')
  in_mat <- rbind(in_mat, df)
}

# run the panConusig
SC_mat <- panConusig::getMatrix(in_mat)

# save the sample-by-component matrix
## output directory
out_dir <- snakemake@params[['outdir']]
output_prefix <- paste0(out_dir, 'panConusig')
write.table(SC_mat, file = paste0(output_prefix, '.SCmatrix.txt'))
saveRDS(SC_mat, file = paste0(output_prefix, '.SCmatrix.rds'), compress = FALSE)

# apply cosine similarity to find the closest signature for each sample
## load the reference signature-by-component matrix
ref_input <- snakemake@params[['def_SC']]
ref_mat <- read.csv(ref_input, sep = '\t', row.names = 1)
panConusig_sig <- colnames(ref_mat)

## normalization on the SC matrix
samples <- colnames(SC_mat)
SC_sum <- colSums(SC_mat)
for (sample_ID in samples) {
  for (i in 1:48) {
    SC_mat[i, sample_ID] <- SC_mat[i, sample_ID] / SC_sum[sample_ID]
  }
}

## calculation
cos_mat <- cbind(ref_mat, SC_mat)
cos_mat <- as.matrix(cos_mat)
cos_mat <- cosine(cos_mat)

# extract the required sample-by-signature similarity matrix
SSmatrix <- cos_mat[samples,]
SSmatrix <- as.data.frame(SSmatrix) %>% select(all_of(panConusig_sig))

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
SSmatrix <- SSmatrix %>% select('sample', 'enrich', 'CN1':'CN25')
row.names(SSmatrix) <- 1:nrow(SSmatrix)
write.table(SSmatrix, file = paste0(output_prefix,'.SSmatrix.tsv'), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
saveRDS(SSmatrix, file = paste0(output_prefix, '.SSmatrix.rds'), compress = FALSE)