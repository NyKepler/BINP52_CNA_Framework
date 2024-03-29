# Title: panConusig_local.R
# Author: Guyuan TANG
# Date: 2024/03/21

# Description: the step 3 in the panConusig_pair_local.R


###### 3. panConusig ######
library(panConusig)
library(tidyverse)
library(lsa)
# load the sample list
sample_df <- read.csv('/home/researcher/TangGY/BINP52/Workflow/Draft/config/sample_pair_panConusig.tsv', sep = '\t')
# prepare the input matrix
in_mat <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(in_mat) <- c('sample', 'chr', 'startpos', 'endpos', 'nMajor', 'nMinor')
for (sampleID in sample_df$Sample) {
  cna_profile <- paste0("/home/researcher/TangGY/BINP52/Workflow/Draft/results/",sampleID,"/06_panConusig/ASCAT_out/",sampleID,"_as_cna_profile.tsv")
  # cna_profile <- paste0('panConusig/',sampleID,'_as_cna_profile.tsv')
  df <- read.csv(cna_profile, sep = '\t')
  in_mat <- rbind(in_mat, df)
}
## remove the NA values in the profiles
in_mat <- filter(in_mat, !is.na(nMajor))

# run the panConusig
SC_mat <- panConusig::getMatrix(in_mat)

# save the sample-by-component matrix
## output directory
out_dir <- "/home/researcher/TangGY/BINP52/Workflow/Draft/results/signatures/panConusig/"
dir.create(out_dir)
# out_dir <- "panConusig/"
output_prefix <- paste0(out_dir, 'panConusig')
saveRDS(SC_mat, file = paste0(output_prefix, '.SCmatrix.rds'), compress = FALSE)
SC_mat_out <- as.data.frame(SC_mat)
SC_mat_out$component <- rownames(SC_mat)
SC_mat_out <- SC_mat_out %>% select('component',colnames(SC_mat))
write.table(SC_mat_out, file = paste0(output_prefix, '.SCmatrix.txt'), sep = '\t', quote = FALSE)

# apply cosine similarity to find the closest signature for each sample
## load the reference signature-by-component matrix
ref_input <- '/home/researcher/TangGY/BINP52/Workflow/Draft/resources/Panconusig_id.txt'
# ref_input <- 'panConusig/Panconusig_id.txt'
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

print('Step 3. panConusig finished.')