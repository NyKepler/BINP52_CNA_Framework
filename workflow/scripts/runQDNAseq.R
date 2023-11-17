#!/usr/bin/env Rscript

# Title: runQDNAseq.R
# Author: Guyuan TANG
# Date: 2023/11/16

# Description: this script will use QDNAseq package to generate relative copy number profile and copy number plot for the input sample.

# Steps:
## 1. Bin annotation
## 2. Processing sorted.dedup.bam file
## 3. Generating relative copy number profile


####### 1. Bin annotation #######
library(QDNAseq)
library(QDNAseq.hg19)

sample_name <- snakemake@params[['sample']]

outdir <- snakemake@params[['outdir']]

# set the fixed bin size
bins <- getBinAnnotations(binSize = snakemake@params[['binsize']])

# read the BAM file and derive raw counts per bin
readCounts <- binReadCounts(bins, bamfiles = snakemake@input[['bamfile']])


####### 2. Processing BAM files #######
# apply filters (blacklist) to the raw read counts per bin
readCounts_Filtered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

# estimate the correction for GC content and mappability
readCounts_Filtered <- estimateCorrection(readCounts_Filtered)
## plot the noise
pdf(paste0(outdir, sample_name, '.sd.pdf'))
noisePlot(readCounts_Filtered)
dev.off()

# run the correction
copyNumbers <- correctBins(readCounts_Filtered)
# run the normalization
copyNumbers_Normalized <- normalizeBins(copyNumbers)
# run the LOESS smoothening
copyNumbers_Smooth <- smoothOutlierBins(copyNumbers_Normalized)

# output the final bins
exportBins(copyNumbers_Smooth, file = snakemake@output[['igv']], format = 'igv', logTransform = FALSE)


####### 3. Generating segmented relative copy number profile #######
# perform segmentation on bins
copyNumbers_Segmented <- segmentBins(copyNumbers_Smooth, transformFun = "sqrt")
# further normalize the segmented bins
copyNumbers_Segmented <- normalizeSegmentedBins(copyNumbers_Segmented)
# plot the segmented counts per bin
pdf(paste0(outdir, sample_name, '_CN_seg.pdf'))
plot(copyNumbers_Segmented)
dev.off()

# export the final segmented bins
exportBins(copyNumbers_Segmented, file = snakemake@output[['seg_tsv']], format = 'tsv', type = 'segments', logTransform = FALSE)
saveRDS(copyNumbers_Segmented, file = snakemake@output[['rds']], compress = FALSE)










