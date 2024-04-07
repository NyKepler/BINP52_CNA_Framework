#!/usr/bin/env Rscript

# Title: get_solutions_samptab.R
# Author: Guyuan TANG
# Date: 2023-12-19

# Description: the script is designed for generating the sample tsv file for second part of the Snakemake workflow in the project. Output tsv file may contain the following columns: sample, patient, type, group, binsize, ploidy, and cellularity.

library(tidyverse)
library(readxl)


# The path of solution xlsx
solution_xlsx <- 'Data/MaNiLa_All_samplesheet_solutions.xlsx'
# The name of the method
method_name <- 'rascal'


# Prepare the final output table
sample_table <- read_excel(solution_xlsx, sheet = 'original')
sample_table <- sample_table[,-grep("Fastq|Sample", colnames(sample_table))]
sample_table <- rename(sample_table, c('Sample' = 'Library')) %>% select(Sample, Patient, Type, Group)
sample_table['Binsize'] <- NA
sample_table['rds'] <- NA
sample_table['Ploidy'] <- NA
sample_table['Cellularity'] <- NA


# Add the bin size and solution for each sample
sample_list <- sample_table$Sample
## add the bin size according to the groups
sample_table[which(sample_table$Group %in% c("ffTumor", "ArchivalVS", "MaNiLaVS")), 'Binsize'] <- 30
sample_table[which(sample_table$Group == "ffpe"), 'Binsize'] <- 100
sample_table$Binsize[is.na(sample_table$Binsize)] <- 50
## use the binsize to extract ploidy and cellularity from the correct sheet
for (ID in sample_list) {
  bins <- sample_table[which(sample_table$Sample == ID), 'Binsize']
  ### QDNAseq RDS file path
  sample_table[which(sample_table$Sample == ID), 'rds'] <- paste0('results/', ID, '/04_relative_CN/', bins, 'kb/', ID, '_', bins, 'kb.rds')
  solution_df <- read_excel(solution_xlsx, sheet = paste0(method_name, '_', bins, 'kb'))
  ### ploidy
  sample_table[which(sample_table$Sample == ID), 'Ploidy'] <- solution_df[which(solution_df$Library == ID), paste0(method_name, '_', bins, 'kb_ploidy_1')]
  ### cellularity
  sample_table[which(sample_table$Sample == ID), 'Cellularity'] <- solution_df[which(solution_df$Library == ID), paste0(method_name, '_', bins, 'kb_cellularity_1')]
}
## remove those without solutions
### Note: 2 samples in ArchivalVS group will be removed
sample_table <- filter(sample_table, Ploidy!=-1 & Cellularity!=-1 )

# change the cellularity 1 into 0
sample_table$Cellularity[which(sample_table$Cellularity == 1)] <- 0

# output the table into tsv file
write.table(sample_table, file = 'config/solution_sample.tsv', sep = '\t', col.names = TRUE, quote = FALSE, row.names = FALSE)




