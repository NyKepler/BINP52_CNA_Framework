# Title: usr_battenberg_pair.R
# Author: Guyuan TANG
# Date: 2024/03/19

# Description: we altered parameters (read depth filter and chromosome names) to fit our sWGS sample data. And the script only includes the first steps for original Battenberg, because we only require those outputs for running further steps in ASCAT.SC.

# Steps:
## 1. prepare_wgs
## 2. run_haplotyping

######### Battenberg main function (new) #########
library(Battenberg)
library(foreach)
library(doParallel)
library(parallel)

battenberg_new = function(analysis="paired", tumourname, normalname, tumour_data_file, normal_data_file, imputeinfofile, g1000prefix, problemloci, gccorrectprefix=NULL,
                      repliccorrectprefix=NULL, g1000allelesprefix=NA, ismale=NA, data_type="wgs", impute_exe="impute2", allelecounter_exe="alleleCounter", nthreads=8, platform_gamma=1, phasing_gamma=1,
                      segmentation_gamma=10, segmentation_kmin=3, phasing_kmin=1, clonality_dist_metric=0, ascat_dist_metric=1, min_ploidy=1.6,
                      max_ploidy=4.8, min_rho=0.1, min_goodness=0.63, uninformative_BAF_threshold=0.51, min_normal_depth=10, min_base_qual=20,
                      min_map_qual=35, calc_seg_baf_option=3, skip_allele_counting=F, skip_preprocessing=F, skip_phasing=F, externalhaplotypefile = NA,
                      usebeagle=FALSE,
                      beaglejar=NA,
                      beagleref.template=NA,
                      beagleplink.template=NA,
                      beaglemaxmem=10,
                      beaglenthreads=1,
                      beaglewindow=40,
                      beagleoverlap=4,
                      javajre="java",
                      write_battenberg_phasing = T, multisample_relative_weight_balanced = 0.25, multisample_maxlag = 100, segmentation_gamma_multisample = 5,
                      snp6_reference_info_file=NA, apt.probeset.genotype.exe="apt-probeset-genotype", apt.probeset.summarize.exe="apt-probeset-summarize",
                      norm.geno.clust.exe="normalize_affy_geno_cluster.pl", birdseed_report_file="birdseed.report.txt", heterozygousFilter="none",
                      prior_breakpoints_file=NULL, GENOMEBUILD="hg19", chrom_coord_file=NULL) {
  
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")
  
  if (analysis == "cell_line"){
    calc_seg_baf_option=1
    phasing_gamma=1
    phasing_kmin=2
    segmentation_gamma=20
    segmentation_kmin=3
    # no matched normal required, but we  are generating normal counts which have this name coded
    normalname = paste0(tumourname, "_normal")
    # other cell_line specific parameter values
  }
  if (analysis == "germline"){
    calc_seg_baf_option=1
    phasing_gamma=3
    phasing_kmin=1
    segmentation_gamma=3
    segmentation_kmin=3
  }
  
  if (data_type=="wgs" & is.na(ismale)) {
    stop("Please provide a boolean denominator whether this sample represents a male donor")
  }
  
  if (data_type=="wgs" & is.na(g1000allelesprefix)) {
    stop("Please provide a path to 1000 Genomes allele reference files")
  }
  
  if (data_type=="wgs" & is.null(gccorrectprefix)) {
    stop("Please provide a path to GC content reference files")
  }
  
  if (!file.exists(problemloci)) {
    stop("Please provide a path to a problematic loci file")
  }
  
  if (!file.exists(imputeinfofile)) {
    stop("Please provide a path to an impute info file")
  }
  
  # check whether the impute_info.txt file contains correct paths
  check.imputeinfofile(imputeinfofile = imputeinfofile, is.male = ismale, usebeagle = usebeagle)
  
  # check whether multisample case
  nsamples <- length(tumourname)
  if (nsamples > 1) {
    if (length(skip_allele_counting) < nsamples) {
      skip_allele_counting = rep(skip_allele_counting[1], nsamples)
    }
    if (length(skip_preprocessing) < nsamples) {
      skip_preprocessing = rep(skip_preprocessing[1], nsamples)
    }
    if (length(skip_phasing) < nsamples) {
      skip_phasing = rep(skip_phasing[1], nsamples)
    }
  }
  
  if (data_type=="wgs" | data_type=="WGS") {
    if (nsamples > 1) {
      print(paste0("Running Battenberg in multisample mode on ", nsamples, " samples: ", paste0(tumourname, collapse = ", ")))
    }
    chrom_names = get.chrom.names.new(imputeinfofile, ismale, analysis=analysis)
  } else if (data_type=="snp6" | data_type=="SNP6") {
    if (nsamples > 1) {
      stop(paste0("Battenberg multisample mode has not been tested with SNP6 data"))
    }
    chrom_names = get.chrom.names(imputeinfofile, TRUE)
    logr_file = paste(tumourname, "_mutantLogR.tab", sep="")
    allelecounts_file = NULL
  }
  print(chrom_names) 
  for (sampleidx in 1:nsamples) {
    
    if (!skip_preprocessing[sampleidx]) {
      if (data_type=="wgs" | data_type=="WGS") {
        # Setup for parallel computing
        clp = parallel::makeCluster(nthreads)
        doParallel::registerDoParallel(clp)
        
        if (analysis == "paired"){
          prepare_wgs_new(chrom_names=chrom_names,
                      tumourbam=tumour_data_file[sampleidx],
                      normalbam=normal_data_file,
                      tumourname=tumourname[sampleidx],
                      normalname=normalname,
                      g1000allelesprefix=g1000allelesprefix,
                      g1000prefix=g1000prefix,
                      gccorrectprefix=gccorrectprefix,
                      repliccorrectprefix=repliccorrectprefix,
                      min_base_qual=min_base_qual,
                      min_map_qual=min_map_qual,
                      allelecounter_exe=allelecounter_exe,
                      min_normal_depth=min_normal_depth,
                      nthreads=nthreads,
                      skip_allele_counting=skip_allele_counting[sampleidx],
                      skip_allele_counting_normal = (sampleidx > 1))
          
        } else if (analysis == "cell_line") {
          Battenberg::prepare_wgs_cell_line(chrom_names=chrom_names,
                                chrom_coord=chrom_coord_file,
                                tumourbam=tumour_data_file,
                                tumourname=tumourname,
                                g1000lociprefix=g1000prefix,
                                g1000allelesprefix=g1000allelesprefix, 
                                gamma_ivd=1e5,
                                kmin_ivd=50,
                                centromere_noise_seg_size=1e6,
                                centromere_dist=5e5,
                                min_het_dist=1e5, 
                                gamma_logr=100,
                                length_adjacent=5e4,
                                gccorrectprefix=gccorrectprefix, 
                                repliccorrectprefix=repliccorrectprefix,
                                min_base_qual=min_base_qual,
                                min_map_qual=min_map_qual, 
                                allelecounter_exe=allelecounter_exe,
                                min_normal_depth=min_normal_depth,
                                skip_allele_counting=skip_allele_counting[sampleidx])
        } else if (analysis == "germline"){
          #prepare_wgs_germline(chrom_names=chrom_names,
          #                      chrom_coord=chrom_coord,
          #                      germlinebam=GERMLINEBAM,
          #                      germlinename=GERMLINENAME,
          #                      g1000lociprefix=G1000PREFIX_AC,
          #                      g1000allelesprefix=G1000PREFIX,
          #                      gamma_ivd=1e5,
          #                      kmin_ivd=50,
          #                      centromere_noise_seg_size=1e6,
          #                      centromere_dist=5e5,
          #                      min_het_dist=2e3,
          #                      gamma_logr=100,
          #                      length_adjacent=5e4,
          #                      gccorrectprefix=GCCORRECTPREFIX,
          #                      repliccorrectprefix=RTCORRECTPREFIX,
          #                      min_base_qual=MIN_BASE_QUAL,
          #                      min_map_qual=MIN_MAP_QUAL,
          #                      allelecounter_exe=ALLELECOUNTER,
          #                      min_normal_depth=MIN_NORMAL_DEPTH,
          #                      skip_allele_counting=F)
        }
        
        
        # Kill the threads
        parallel::stopCluster(clp)
        
      } else if (data_type=="snp6" | data_type=="SNP6") {
        
        Battenberg::prepare_snp6(tumour_cel_file=tumour_data_file[sampleidx],
                     normal_cel_file=normal_data_file,
                     tumourname=tumourname[sampleidx],
                     chrom_names=chrom_names,
                     snp6_reference_info_file=snp6_reference_info_file,
                     apt.probeset.genotype.exe=apt.probeset.genotype.exe,
                     apt.probeset.summarize.exe=apt.probeset.summarize.exe,
                     norm.geno.clust.exe=norm.geno.clust.exe,
                     birdseed_report_file=birdseed_report_file)
        
      } else {
        print("Unknown data type provided, please provide wgs or snp6")
        q(save="no", status=1)
      }
    }
    
    if (data_type=="snp6" | data_type=="SNP6") {
      # Infer what the gender is - WGS requires it to be specified
      gender = infer_gender_birdseed(birdseed_report_file)
      ismale = gender == "male"
    }
    
    
    if (!skip_phasing[sampleidx]) {
      
      # if external phasing data is provided (as a vcf), split into chromosomes for use in haplotype reconstruction
      if (!is.na(externalhaplotypefile) && file.exists(externalhaplotypefile)) {
        externalhaplotypeprefix <- paste0(normalname, "_external_haplotypes_chr")
        
        # if these files exist already, no need to split again
        if (any(!file.exists(paste0(externalhaplotypeprefix, 1:length(chrom_names), ".vcf")))) {
          
          print(paste0("Splitting external phasing data from ", externalhaplotypefile))
          Battenberg::split_input_haplotypes(chrom_names = chrom_names,
                                 externalhaplotypefile = externalhaplotypefile,
                                 outprefix = externalhaplotypeprefix)
        } else {
          print("No need to split, external haplotype files per chromosome found")
        }
      } else {
        externalhaplotypeprefix <- NA
      }
      
      # Setup for parallel computing
      clp = parallel::makeCluster(nthreads)
      doParallel::registerDoParallel(clp)
      
      # Reconstruct haplotypes
      # mclapply(1:length(chrom_names), function(chrom) {
      foreach::foreach (i=1:length(chrom_names)) %dopar% {
        chrom = chrom_names[i]
        print(chrom)
        
        Battenberg::run_haplotyping(chrom=chrom,
                        tumourname=tumourname[sampleidx],
                        normalname=normalname,
                        ismale=ismale,
                        imputeinfofile=imputeinfofile,
                        problemloci=problemloci,
                        impute_exe=impute_exe,
                        min_normal_depth=min_normal_depth,
                        chrom_names=chrom_names,
                        snp6_reference_info_file=snp6_reference_info_file,
                        heterozygousFilter=heterozygousFilter,
                        usebeagle=usebeagle,
                        beaglejar=beaglejar,
                        beagleref=gsub("CHROMNAME", chrom, beagleref.template),
                        beagleplink=gsub("CHROMNAME", chrom, beagleplink.template),
                        beaglemaxmem=beaglemaxmem,
                        beaglenthreads=beaglenthreads,
                        beaglewindow=beaglewindow,
                        beagleoverlap=beagleoverlap,
                        externalhaplotypeprefix=externalhaplotypeprefix,
                        use_previous_imputation=(sampleidx > 1))
      }#, mc.cores=nthreads)
      
      # Kill the threads as from here its all single core
      parallel::stopCluster(clp) }
  }
}


## the new prepare_wgs_new function
prepare_wgs_new = function(chrom_names, tumourbam, normalbam, tumourname, normalname, g1000allelesprefix, g1000prefix, gccorrectprefix,
                       repliccorrectprefix, min_base_qual, min_map_qual, allelecounter_exe, min_normal_depth, nthreads, skip_allele_counting, skip_allele_counting_normal = F) {
  
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")
  
  if (!skip_allele_counting) {
    # Obtain allele counts for 1000 Genomes locations for both tumour and normal
    foreach::foreach(i=1:length(chrom_names)) %dopar% {
      Battenberg::getAlleleCounts(bam.file=tumourbam,
                      output.file=paste(tumourname,"_alleleFrequencies_chr", chrom_names[i], ".txt", sep=""),
                      g1000.loci=paste(g1000prefix, chrom_names[i], ".txt", sep=""),
                      min.base.qual=min_base_qual,
                      min.map.qual=min_map_qual,
                      allelecounter.exe=allelecounter_exe)
      
      if (!skip_allele_counting_normal) {
        Battenberg::getAlleleCounts(bam.file=normalbam,
                        output.file=paste(normalname,"_alleleFrequencies_chr", chrom_names[i], ".txt",  sep=""),
                        g1000.loci=paste(g1000prefix, chrom_names[i], ".txt", sep=""),
                        min.base.qual=min_base_qual,
                        min.map.qual=min_map_qual,
                        allelecounter.exe=allelecounter_exe)
      }
    }
  }
  
  # remove the 'chr' in the allele frequency output files
  Battenberg::standardiseChrNotation(tumourname,normalname)
  
  # Obtain BAF and LogR from the raw allele counts
  Battenberg::getBAFsAndLogRs(tumourAlleleCountsFile.prefix=paste(tumourname,"_alleleFrequencies_chr", sep=""),
                  normalAlleleCountsFile.prefix=paste(normalname,"_alleleFrequencies_chr", sep=""),
                  figuresFile.prefix=paste(tumourname, "_", sep=''),
                  BAFnormalFile=paste(tumourname,"_normalBAF.tab", sep=""),
                  BAFmutantFile=paste(tumourname,"_mutantBAF.tab", sep=""),
                  logRnormalFile=paste(tumourname,"_normalLogR.tab", sep=""),
                  logRmutantFile=paste(tumourname,"_mutantLogR.tab", sep=""),
                  combinedAlleleCountsFile=paste(tumourname,"_alleleCounts.tab", sep=""),
                  chr_names=chrom_names,
                  g1000file.prefix=g1000allelesprefix,
                  minCounts=min_normal_depth,
                  samplename=tumourname)
  # Perform GC correction
  Battenberg::gc.correct.wgs(Tumour_LogR_file=paste(tumourname,"_mutantLogR.tab", sep=""),
                 outfile=paste(tumourname,"_mutantLogR_gcCorrected.tab", sep=""),
                 correlations_outfile=paste(tumourname, "_GCwindowCorrelations.txt", sep=""),
                 gc_content_file_prefix=gccorrectprefix,
                 replic_timing_file_prefix=repliccorrectprefix,
                 chrom_names=chrom_names)
}





## required check.imputeinfofile
check.imputeinfofile = function(imputeinfofile, is.male, usebeagle) {
  impute.info = parse.imputeinfofile(imputeinfofile, is.male)
  if (usebeagle){
    if (any(!file.exists(impute.info$impute_legend))) {
      print("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
      stop("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
    }
  } else {
    if (any(!file.exists(impute.info$impute_legend) | !file.exists(impute.info$genetic_map) | !file.exists(impute.info$impute_hap))) {
      print("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
      stop("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
    }
  }
}

# new get.chrom.names function to excclude X 
get.chrom.names.new = function(imputeinfofile, is.male, chrom=NA, analysis="paired") {
  chrom_names = unique(parse.imputeinfofile(imputeinfofile, is.male, chrom=chrom)$chrom)
  if (analysis=="cell_line" | analysis=="germline" | analysis=="paired") {
    # Both cell line and germline analysis do not yield usable data on X and Y, so remove
    chrom_names = chrom_names[!chrom_names %in% c("X", "Y")]
  }
  return(chrom_names)
}
