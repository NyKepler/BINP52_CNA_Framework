#!/usr/bin/env Rscript

# Title: fit_absoluteCN.R
# Author: Guyuan TANG
# Date: 2023/12/27

# Description: this script will be used for generating the absolute copy number profile (including copy numbers and segments) for each sample. It is designed based on the fit_absolute_copy_number.R provided by the rascal package.

# Steps:
## 1. take in the sample table and extract the information (binsize, ploidy and cellularity) of the selected sample;
## 2. load the corresponding QDNAseq RDS file and use the information to calculate the absolute copy numbers, and apply segmentation on the absolute copy number profile;
## 3. export the absolute copy numbers and segments for the following analysis in the workflow.

library(tibble)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(QDNAseq)
library(rascal)

# specify the output location
outdir <- snakemake@params[['outdir']]

####### 1. Sample Information #######
# load the sample table
samptable_path <- snakemake@input[['sample_table']]
sample_df <- read.table(samptable_path, sep = '\t', header = 1)
# the information of the selected sample
sample_ID <- snakemake@params[['sample']]
bins <- sample_df[which(sample_df$Sample == sample_ID), 'Binsize']
ploidy <- as.numeric(sample_df[which(sample_df$Sample == sample_ID), 'Ploidy'])
cellularity <- as.numeric(sample_df[which(sample_df$Sample == sample_ID), 'Cellularity'])
rds <- sample_df[which(sample_df$Sample == sample_ID), 'rds']

####### 2. Calculate Absolute Copy Number #######
# function for extracting the copy number for a given sample
# can handle QDNAseqCopyNumbers object or a copy number data frame
copy_number_for_sample <- function(copy_number, sample) {
  if (any(class(copy_number) == "QDNAseqCopyNumbers")) {
    copy_number <- copy_number[,sample]
    copy_number_values <- Biobase::assayDataElement(copy_number, "copynumber")[,1]
    segmented_values <- Biobase::assayDataElement(copy_number, "segmented")[,1]
    Biobase::fData(copy_number) %>%
      rownames_to_column(var = "id") %>%
      as_tibble() %>%
      select(id, chromosome, start, end) %>%
      mutate(across(c(start, end), as.integer)) %>%
      mutate(chromosome = factor(chromosome, levels = unique(chromosome))) %>%
      mutate(sample = sample) %>%
      mutate(copy_number = copy_number_values) %>%
      mutate(segmented = segmented_values) %>%
      select(sample, chromosome, start, end, copy_number, segmented)
  } else {
    filter(copy_number, sample == !!sample)
  }
}

# load the QDNAseq rds file
copy_number <- readRDS(rds)
# check contents are correct and obtain sample names
if (any(class(copy_number) == "data.frame")) {
  required_columns <- c("sample", "chromosome", "start", "end", "segmented")
  missing_columns <- setdiff(required_columns, colnames(copy_number))
  if (length(missing_columns) > 0) stop("missing columns in", input_file, ": ", str_c(missing_columns, collapse = ", "))
  samples <- sort(unique(copy_number$sample))
} else if (any(class(copy_number) == "QDNAseqCopyNumbers")) {
  if (!requireNamespace("QDNAseq", quietly = TRUE)) stop("QDNAseq package needs to be installed")
  samples <- sort(Biobase::sampleNames(copy_number))
} else {
  stop(rds, " should contain either a data frame or a QDNAseqCopyNumbers object")
}

# absolute copy number calculation
sample_copy_number <- copy_number_for_sample(copy_number, sample = paste0(sample_ID, '.sorted.dedup'))
relative_copy_number <- mutate(sample_copy_number, across(c(copy_number, segmented), ~ . / median(segmented, na.rm = TRUE)))
segments <- copy_number_segments(relative_copy_number)
absolute_copy_number <- mutate(relative_copy_number, 
                               across(c(copy_number, segmented),
                                      relative_to_absolute_copy_number, ploidy, cellularity))
absolute_segments <- copy_number_segments(absolute_copy_number)
# rename the column and add the sample column
absolute_copy_number <- rename(absolute_copy_number, c('segVal'='segmented'))
absolute_copy_number['sample'] <- sample_ID
absolute_segments <- rename(absolute_segments, c('segVal'='copy_number'))
absolute_segments['sample'] <- sample_ID

####### 3. Export Absolute Copy Number #######
output_CN <- paste0(outdir, sample_ID, '_', bins, 'kb_CN.tsv')
write.table(absolute_copy_number, file = output_CN, row.names = FALSE, sep = '\t', quote = FALSE)
output_seg <- paste0(outdir, sample_ID, '_', bins, 'kb_seg.tsv')
write.table(absolute_segments, file = output_seg, row.names = FALSE, sep = '\t', quote = FALSE)







