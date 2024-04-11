# Title: panConusig_pair_local.R
# Author: Guyuan TANG
# Date: 2024/03/21

# Description: to avoid errors, we ran this local R scripts on compuster instead of running it within Snakemake workflow.

# Steps:
## 1. run Battenberg to generate the allele frequency files and the phased files
## 2. run ASCAT.sc to generate the allele-specific copy number profiles
## 3. run panConusig to generate the sample-by-component matrix and the sample-by-signatures (cosine similarity) matrix

library(ASCAT.sc)
library(tidyverse)

###### 2. run ASCAT.sc ######
# load the sample table
sample_df <- read.csv('/home/researcher/TangGY/BINP52/Workflow/Draft/config/sample_pair_panConusig.tsv', sep = '\t')

for (sampleID in sample_df$Sample) {
  print(sampleID)
  # set the working directory
  working_dir <- paste0("/home/researcher/TangGY/BINP52/Workflow/Draft/results/",sampleID,"/06_panConusig/")
  setwd(working_dir)
  
  # run ASCAT.sc
  ASCAT_out <- paste0("/home/researcher/TangGY/BINP52/Workflow/Draft/results/",sampleID,"/06_panConusig/ASCAT_out/")
  dir.create(ASCAT_out)
  results_output = paste0("/home/researcher/TangGY/BINP52/Workflow/Draft/results/",sampleID,'/06_panConusig/')
  ## path to required files
  path_to_phases=lapply(sampleID, function(SAMPLENAME)
  {
    paste0(results_output, SAMPLENAME,
           "_beagle5_output_chr",c(1:22),".txt.vcf.gz")
  })
  
  list_ac_counts_paths=lapply(sampleID, function(SAMPLENAME)
  {
    paste0(results_output, SAMPLENAME,
           "_alleleFrequencies_chr",c(1:22),".txt")
  })
  bins <- sample_df[which(sample_df$Sample==sampleID),'Binsize']
  bins <- as.numeric(bins) * 1000
  TUMOURBAM = paste0("/home/researcher/TangGY/BINP52/Workflow/Draft/results/",sampleID,'/03_clean_up/',sampleID,'.sorted.dedup.bam')
  ## run the main programme
  res <- run_sc_sequencing(tumour_bams=TUMOURBAM,
                           allchr=paste0("chr",c(1:22)),
                           sex='female',
                           binsize=bins,
                           chrstring_bam="chr",
                           purs = seq(0.01, 1, 0.001), #if start from 0, will affect the following operations
                           ploidies = seq(1.7, 5, 0.01),
                           maxtumourpsi=5,
                           build="hg19",
                           MC.CORES=30,
                           outdir=ASCAT_out,
                           projectname="ASCAT_CN",
                           segmentation_alpha=0.01,
                           path_to_phases=path_to_phases,
                           list_ac_counts_paths=list_ac_counts_paths,
                           predict_refit=TRUE,
                           multipcf=FALSE)
  
  ## alter the format of the final output
  in_file <- paste0(ASCAT_out, 'as_cna_profile_', sampleID, '.sorted.dedup.bam_bam1.txt')
  df <- read.csv(in_file, sep = '\t')
  df$sample <- sampleID
  df <- df %>% select('sample', 'chr', 'startpos', 'endpos', 'nA', 'nB')
  colnames(df) <- c('sample', 'chr', 'startpos', 'endpos', 'nMajor', 'nMinor')
  # print the file
  write.table(df, file = paste0(ASCAT_out, sampleID, '_as_cna_profile.tsv'), quote = FALSE, sep = '\t', row.names = FALSE)
}
print('Step 2. ASCAT.sc finished.')