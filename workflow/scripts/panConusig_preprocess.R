# Title: panConusig_preprocess.R
# Author: Guyuan TANG
# Date: 2024-03-06

# Description: this is used for the preprocessing steps before generating the panConusig. We used Battenberg to perform the haplotyping and ASCAT.sc to generate the allele-specific copy number profiles.

# Steps:
## 1. source the 'usr_battenberg.R' to alter the depth filter to 2 (the original was 10)
## 2. run new Battenberg function to generate allele frequency files and phased files; remove the files that are not required for next step
## 3. run ASCAT.sc on the allele frequency files and phased files
## 4. alter the format of the output files from ASCAT.sc to meet the requirement for the next step in Snakemake workflow

library(Battenberg)
library(ASCAT.sc)

##### 1. user version of battenberg #####
usr_battenberg <- snakemake@input[['usr_battenberg']]
source(usr_battenberg)

setwd(snakemake@params[['workdir']])

##### 2. run Battenberg #####
sampleID <- snakemake@params[['sampleID']]
TUMOURNAME <- sampleID
TUMOURBAM <- snakemake@input[['bam_file']]
IS.MALE = FALSE
SKIP_ALLELECOUNTING = FALSE
SKIP_PREPROCESSING = FALSE
SKIP_PHASING = FALSE
NTHREADS = snakemake@threads
PRIOR_BREAKPOINTS_FILE = NULL

analysis <- 'cell_line'
GENOMEBUILD <- "hg19"
USEBEAGLE <- TRUE
JAVAJRE <- "java"
ALLELECOUNTER <- "alleleCounter"
IMPUTE_EXE <- "impute2"

# location of the reference files
if (GENOMEBUILD=="hg19") {
  impute_basedir = snakemake@params[['impute_ref_dir']]
  IMPUTEINFOFILE = file.path(impute_basedir, "battenberg_impute_v3/impute_info.txt")
  G1000ALLELESPREFIX = file.path(impute_basedir, "battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr")
  G1000LOCIPREFIX = file.path(impute_basedir, "battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr")
  GCCORRECTPREFIX = file.path(impute_basedir, "battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_")
  REPLICCORRECTPREFIX = file.path(impute_basedir, "battenberg_wgs_replic_correction_1000g_v3/1000_genomes_replication_timing_chr_")
  
  # WGS specific static
  PROBLEMLOCI = file.path(impute_basedir, "probloci_270415.txt.gz")
  GENOME_VERSION = "b37"
  GENOMEBUILD = "hg19"
  BEAGLE_BASEDIR = snakemake@params[['beagle_ref_dir']]
  BEAGLEJAR = file.path(BEAGLE_BASEDIR, "beagle.22Jul22.46e.jar")
  BEAGLEREF.template = file.path(BEAGLE_BASEDIR, GENOME_VERSION, "chrCHROMNAME.1kg.phase3.v5a.b37.bref3")
  BEAGLEPLINK.template = file.path(BEAGLE_BASEDIR, GENOME_VERSION, "plink.chrCHROMNAME.GRCh37.map")
  
  CHROM_COORD_FILE = file.path(impute_basedir, "gcCorrect_chromosome_coordinates_hg19.txt")
  
  PLATFORM_GAMMA = 1
  PHASING_GAMMA = 1
  SEGMENTATION_GAMMA = 10
  SEGMENTATIIN_KMIN = 3
  PHASING_KMIN = 1
  CLONALITY_DIST_METRIC = 0
  ASCAT_DIST_METRIC = 1
  MIN_PLOIDY = 1
  MAX_PLOIDY = 4.8
  MIN_RHO = 0
  MIN_GOODNESS_OF_FIT = 0.63
  BALANCED_THRESHOLD = 0.51
  MIN_NORMAL_DEPTH = 10
  MIN_BASE_QUAL = 20
  MIN_MAP_QUAL = 30
  CALC_SEG_BAF_OPTION = 1
  
}

battenberg_new(analysis=analysis,
               tumourname=TUMOURNAME, 
               normalname=NA, 
               tumour_data_file=TUMOURBAM, 
               normal_data_file=NA, 
               ismale=IS.MALE, 
               imputeinfofile=IMPUTEINFOFILE, 
               g1000prefix=G1000LOCIPREFIX, 
               g1000allelesprefix=G1000ALLELESPREFIX, 
               gccorrectprefix=GCCORRECTPREFIX, 
               repliccorrectprefix=REPLICCORRECTPREFIX, 
               problemloci=PROBLEMLOCI, 
               data_type="wgs",
               impute_exe=IMPUTE_EXE,
               allelecounter_exe=ALLELECOUNTER,
               usebeagle=USEBEAGLE,
               beaglejar=BEAGLEJAR,
               beagleref=BEAGLEREF.template,
               beagleplink=BEAGLEPLINK.template,
               beaglemaxmem=10,
               beaglenthreads=1,
               beaglewindow=40,
               beagleoverlap=4,
               javajre=JAVAJRE,
               nthreads=NTHREADS,
               platform_gamma=PLATFORM_GAMMA,
               phasing_gamma=PHASING_GAMMA,
               segmentation_gamma=SEGMENTATION_GAMMA,
               segmentation_kmin=SEGMENTATIIN_KMIN,
               phasing_kmin=PHASING_KMIN,
               clonality_dist_metric=CLONALITY_DIST_METRIC,
               ascat_dist_metric=ASCAT_DIST_METRIC,
               min_ploidy=MIN_PLOIDY,
               max_ploidy=MAX_PLOIDY,
               min_rho=MIN_RHO,
               min_goodness=MIN_GOODNESS_OF_FIT,
               uninformative_BAF_threshold=BALANCED_THRESHOLD,
               min_normal_depth=MIN_NORMAL_DEPTH,
               min_base_qual=MIN_BASE_QUAL,
               min_map_qual=MIN_MAP_QUAL,
               calc_seg_baf_option=CALC_SEG_BAF_OPTION,
               skip_allele_counting=SKIP_ALLELECOUNTING,
               skip_preprocessing=SKIP_PREPROCESSING,
               skip_phasing=SKIP_PHASING,
               prior_breakpoints_file=PRIOR_BREAKPOINTS_FILE,
               GENOMEBUILD=GENOMEBUILD,
               chrom_coord_file=CHROM_COORD_FILE)
print('Battenberg analysis finished.')
# prepare the required files in the corresponding sample folder
result_dir <- snakemake@params[['result_dir']]
## copy the required files
for (i in 1:22) {
  allele_freq_file <- paste0(sampleID, "_alleleFrequencies_chr", i, ".txt")
  phased_file <- paste0(sampleID, "_beagle5_output_chr", i, ".txt.vcf.gz")
  file.copy(from = allele_freq_file, to = paste0(result_dir, allele_freq_file))
  file.copy(from = phased_file, to = paste0(result_dir, phased_file))
}
## remove the non-required files (not combine with the copy step to avoid corruptions)
unlink('PCF_plots/', recursive = TRUE)
remove_files <- c(paste0(sampleID, '_mutantBAF.tab'), paste0(sampleID, '_mutantLogR.tab'), paste0(sampleID, '_mutantLogR_gcCorrected.tab'),
                  paste0(sampleID, '_GCwindowCorrelations_beforeCorrection.txt'), paste0(sampleID, '_GCwindowCorrelations_afterCorrection.txt'))
file.remove(remove_files)
for (i in 1:22) {
  chr_remove_file <- c(paste0(sampleID, '_alleleFrequencies_chr', i, '.txt'), paste0(sampleID, '_normal_alleleFrequencies_chr', i, '.txt'),
                       paste0(sampleID, '_beagle5_input_chr', i, '.txt'), paste0(sampleID, '_beagle5_output_chr', i, '.txt.log'), 
                       paste0(sampleID, '_beagle5_output_chr', i, '.txt.vcf.gz'), 
                       paste0(sampleID, '_chr', i, '_heterozygousData.png'), paste0(sampleID, '_chr', i, '_heterozygousMutBAFs_haplotyped.txt'),
                       paste0(sampleID, '_impute_input_chr', i, '.txt'), paste0(sampleID, '_impute_output_chr', i, '_allHaplotypeInfo.txt'))
  file.remove(chr_remove_file)
}


##### 3. run ASCAT.sc on the allele frequency files and phased haplotype files #####
path_to_phases=lapply(sampleID, function(SAMPLENAME)
{
  paste0(result_dir, SAMPLENAME,
         "_beagle5_output_chr",c(1:22),".txt.vcf.gz")
})

list_ac_counts_paths=lapply(sample_ID, function(SAMPLENAME)
{
  paste0(result_dir, SAMPLENAME,
         "_alleleFrequencies_chr",c(1:22),".txt")
})

# output dir
ASCAT_out <- snakemake@params[['ASCAT_outdir']]

# run ASCAT.sc
bins <- as.numeric(snakemake@parmas[['binsize']])
bins <- bins * 1000
res <- run_sc_sequencing(tumour_bams=TUMOURBAM,
                         allchr=paste0("chr",c(1:22)),
                         sex='female',
                         binsize=bins,
                         chrstring_bam="chr",
                         purs = seq(0.01, 1, 0.01), #if start from 0, will affect the following operations
                         ploidies = seq(1,5, 0.01),
                         maxtumourpsi=5,
                         build="hg19",
                         MC.CORES=1,
                         outdir=ASCAT_out,
                         projectname="ASCAT_CN",
                         segmentation_alpha=0.01,
                         path_to_phases=path_to_phases,
                         list_ac_counts_paths=list_ac_counts_paths,
                         predict_refit=TRUE,
                         multipcf=FALSE)


##### 4. alter the format of the final output #####
library(tidyverse)
in_file <- paste0(ASCAT_out, 'as_cna_profile_', sampleID, '.sorted.dedup.bam_bam1.txt')
df <- read.csv(in_file, sep = '\t')
df$sample <- sampleID
df <- df %>% select('sample', 'chr', 'startpos', 'endpos', 'nA', 'nB') %>% rename('nMajor'='nA', 'nMinor'='nB')
# print the file
write.csv(df, file = paste0(ASCAT_out, sampleID, '_as_cna_profile.tsv'), quote = FALSE, sep = '\t')
