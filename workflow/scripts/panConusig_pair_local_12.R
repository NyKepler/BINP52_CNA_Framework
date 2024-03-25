# Title: panConusig_pair_local.R
# Author: Guyuan TANG
# Date: 2024/03/21

# Description: to avoid errors, we ran this local R scripts on compuster instead of running it within Snakemake workflow.

# Steps:
## 1. run Battenberg to generate the allele frequency files and the phased files
## 2. run ASCAT.sc to generate the allele-specific copy number profiles
## 3. run panConusig to generate the sample-by-component matrix and the sample-by-signatures (cosine similarity) matrix

library(Battenberg)
library(ASCAT.sc)
library(tidyverse)


##### 1. preprocess in Battenberg #####
# load the new adjusted function
source("/home/researcher/TangGY/BINP52/Workflow/Draft/workflow/scripts/usr_battenberg_pair.R")
# load the sample table
sample_df <- read.csv('/home/researcher/TangGY/BINP52/Workflow/Draft/config/sample_pair_panConusig.tsv', sep = '\t')

# run battenberg
for (sampleID in sample_df$Sample) {
  # set the working directory
  working_dir <- paste0("/home/researcher/TangGY/BINP52/Workflow/Draft/results/",sampleID,"/06_panConusig/")
  dir.create(working_dir)
  setwd(working_dir)
  # load the parameters
  TUMOURNAME = sampleID
  NORMALNAME = sample_df[which(sample_df$Sample==sampleID),'Sample_ref']
  results_output = "/home/researcher/TangGY/BINP52/Workflow/Draft/results/"
  NORMALBAM = paste0(results_output,NORMALNAME,'/03_clean_up/',NORMALNAME,'.sorted.dedup.bam')
  TUMOURBAM = paste0(results_output,TUMOURNAME,'/03_clean_up/',TUMOURNAME,'.sorted.dedup.bam')
  IS.MALE = FALSE
  SKIP_ALLELECOUNTING = FALSE
  SKIP_PREPROCESSING = FALSE
  SKIP_PHASING = FALSE
  NTHREADS = 30
  PRIOR_BREAKPOINTS_FILE = NULL
  
  analysis = "paired"
  
  JAVAJRE = "java"
  ALLELECOUNTER = "alleleCounter"
  IMPUTE_EXE = "impute2"
  
  GENOMEBUILD = "hg19"
  USEBEAGLE = T
  
  # General static
  # location of the reference files
  if (GENOMEBUILD=="hg19") {
    impute_basedir = "/home/researcher/TangGY/BINP52/Workflow/Draft/resources/battenberg"
    IMPUTEINFOFILE = file.path(impute_basedir, "battenberg_impute_v3/impute_info.txt")
    G1000ALLELESPREFIX = file.path(impute_basedir, "battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr")
    G1000LOCIPREFIX = file.path(impute_basedir, "battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr")
    GCCORRECTPREFIX = file.path(impute_basedir, "battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_")
    REPLICCORRECTPREFIX = file.path(impute_basedir, "battenberg_wgs_replic_correction_1000g_v3/1000_genomes_replication_timing_chr_")
    
    # WGS specific static
    PROBLEMLOCI = file.path(impute_basedir, "probloci_270415.txt.gz")
    GENOME_VERSION = "b37"
    GENOMEBUILD = "hg19"
    BEAGLE_BASEDIR = "/home/researcher/TangGY/BINP52/Workflow/Draft/resources/battenberg/beagle"
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
    MIN_PLOIDY = 1.6
    MAX_PLOIDY = 4.8
    MIN_RHO = 0.1
    MIN_GOODNESS_OF_FIT = 0.63
    BALANCED_THRESHOLD = 0.51
    MIN_NORMAL_DEPTH = 2 #default setting is 10 which is not suitable for sWGS, and for most positions, the good read depth is only 1
    MIN_BASE_QUAL = 20
    MIN_MAP_QUAL = 35
    CALC_SEG_BAF_OPTION = 1
  }
  # run battenberg
  battenberg_new(analysis=analysis,
                 tumourname=TUMOURNAME, 
                 normalname=NORMALNAME, 
                 tumour_data_file=TUMOURBAM, 
                 normal_data_file=NORMALBAM, 
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
                 beaglenthreads=30,
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
  
  # remove the files that are not required
  ## copy the blood samples' allele frequencies files to the correspond directory
  ref_dir <- paste0('/home/researcher/TangGY/BINP52/Workflow/Draft/results/',NORMALNAME,'/06_panConusig/')
  if (!dir.exists(ref_dir)){
    dir.create(ref_dir) }
  for (chr in 1:22) {
    allelefreq_file <- paste0(NORMALNAME, '_alleleFrequencies_chr',chr,'.txt')
    file.copy(from = allelefreq_file, to = paste0(ref_dir,allelefreq_file))
    file.remove(allelefreq_file)
  }
  ## the chromosome files that are not required
  for (chr in 1:22) {
    beagle_input <- paste0(TUMOURNAME,'_beagle5_input_chr',chr,'.txt')
    beagle_out <- paste0(TUMOURNAME,'_beagle5_output_chr',chr,'.txt.log')
    hData <- paste0(TUMOURNAME,'_chr',chr,'_heterozygousData.png')
    hMutBAF <- paste0(TUMOURNAME,'_chr',chr,'_heterozygousMutBAFs_haplotyped.txt')
    impute_input <- paste0(TUMOURNAME,'_impute_input_chr',chr,'.txt')
    impute_output <- paste0(TUMOURNAME,'_impute_output_chr',chr,'_allHaplotypeInfo.txt')
    file.remove(beagle_input,beagle_out,hData,impute_input,impute_output)
  }
  ## other files that are not required
  mutantBAF <- paste0(TUMOURNAME,'_mutantBAF.tab')
  mutantLogR <- paste0(TUMOURNAME,'_mutantLogR.tab')
  mutantLogR_gc <- paste0(TUMOURNAME,'_mutantLogR_gcCorrected.tab')
  normalBAF <- paste0(TUMOURNAME,'_normalBAF.tab')
  normalLogR <- paste0(TUMOURNAME,'_normalLogR.tab')
  GC_before <- paste0(TUMOURNAME,'_GCwindowCorrelations_beforeCorrection.txt')
  GC_after <- paste0(TUMOURNAME,'_GCwindowCorrelations_afterCorrection.txt')
  ac_tab <- paste0(TUMOURNAME,'_alleleCounts.tab')
  file.remove(mutantBAF,mutantLogR,mutantLogR_gc,normalBAF,normalLogR,GC_before, GC_after,ac_tab)

}
print('Step 1. Battenberg on all samples finished.')



###### 2. run ASCAT.sc ######
# load the sample table
sample_df <- read.csv('/home/researcher/TangGY/BINP52/Workflow/Draft/config/sample_pair_panConusig.tsv', sep = '\t')

# remove the intermediate files that are not required
for (sampleID in c('CS2_176')) {
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
                           purs = seq(0.01, 1, 0.01), #if start from 0, will affect the following operations
                           ploidies = seq(1,5, 0.01),
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

###### 3. panConusig ######


