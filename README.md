# README for the CNA Framework
Author: Guyuan Tang  
Date: 2023/10/16 - 

## 1. Description 
### 1.1 The project
This project is designed within a master thesis (BINP52). We aim to develop a pipeline to generate copy number profiles from shallow whole genome sequening (sWGS) samples.  

### 1.2 The pipeline
The pipeline will include steps from preprocessing raw reads with QC to obtaining absolute copy number profiles.  
#### Steps
- (1) Preprocessing. This step includes quality assessment and quality trimming on the raw reads. (`Fastp` will be used for QC and trimming, together with `multiQC` to generate the QC reports.) [Or `fastQC` + `Trimmomatic` + `multiQC`?]
- (2) Alignment. The human reference genome will be indexed. And the reads will be mapped to the reference genome. (`BWA` will be used for both indexing and alignment.)
- (3) Clean-up. After alignment, the SAM files will be sorted and the PCR duplicates will be marked and removed. Also, the .sorted.deduplicated.sam will be converted to BAM files. The BAM files will be indexed for later analysis. (`Picard` will be used for sorting SAM, marking duplicates, removing duplicates and converting SAM to BAM. `samtools` will be used for indexing the BAM files.)
- (4) Relative copy number profile. The BAM files will be analyzed through fixed-size binning, filtering, correction, normalization to generate the read counts per bin. This data will then used for segmentation of bins and generating the relative copy number profile. (`QDNAseq` will be used for this step.)
- (5) Absolute copy number profile. The output file from `QDNAseq` contains relative copy number, and we need to estimate ploidy and cellularity in our samples to generate our final absolute copy number profile for comparison. (`Rascal` will be used for this step.)

## 2. Workflow Details
### 2.1 Preprocessing
This step includes quality assessment and quality trimming on the raw reads.