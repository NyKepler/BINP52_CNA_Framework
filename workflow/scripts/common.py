# Title: common.py
# Author: Guyuan TANG
# Date: 2023/12/29

# Description: the script to store functions used for generating the input and output file lists. Because some file names will have different variables (bin size) which are difficult to only use expand() function in Snakemake.


## 1. generate the file list for output of absolute copy number profiles
def get_output_absolute(sample_df, results):
    files = []
    for sample_name in sample_df.Sample:
        binsize = sample_df.loc[sample_name,'Binsize']
        binsize = binsize + 'kb'
        file_name = results + sample_name + '/05_absolute_CN/' + sample_name + '_' + binsize + '_seg.tsv'
        files.append(file_name)
    return files

## 2. generate the reference file names for Battenberg
def get_ref_battenberg(working_dir):
    files = []
    dir_prefix = working_dir + 'resources/battenberg/'
    for i in range(1,24):
        i = str(i)
        # 1000 genome loci files
        genomelociAllele =  dir_prefix +'battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr' + i + '.txt'
        genomeloci = dir_prefix + 'battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr' + i + '.txt'
        # wgs GC correction files
        wgs_gc_correction = dir_prefix + 'battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_' + i + '.txt'
        # wgs replication correction files
        wgs_rep_correction = dir_prefix + 'battenberg_wgs_replic_correction_1000g_v3/1000_genomes_replication_timing_chr_' + i + '.txt.gz'
        files.extend((genomelociAllele, genomeloci, wgs_gc_correction, wgs_rep_correction))
    impute_list = list(range(1,23))
    impute_list.extend(('X_nonPAR', 'X_PAR1', 'X_PAR2'))
    files.extend((dir_prefix+'battenberg_impute_v3/ALL_1000G_phase1integrated_v3.sample', dir_prefix+'impute_info.txt', dir_prefix+'README_1000G_phase1integrated_v3.txt'))
    for chr in impute_list:
        chr = str(chr)
        # ALL phased haplotypes
        ALL_hap = dir_prefix + 'battenberg_impute_v3/ALL_1000G_phase1integrated_v3_chr' + chr + '_impute.hap.gz'
        ALL_hap_legend = dir_prefix + 'battenberg_impute_v3/ALL_1000G_phase1integrated_v3_chr' + chr + '_impute.legend'
        # genetic map
        genetic_map = dir_prefix + 'battenberg_impute_v3/genetic_map_chr' + chr + '_combined_b37.txt'
        files.extend((ALL_hap, ALL_hap_legend, genetic_map))
    files.extend((dir_prefix+'gcCorrect_chromosome_coordinates_hg19.txt', dir_prefix+'probloci_270415.txt.gz'))
    return files

## 3. generate the reference file names for beagle
def get_ref_beagle(working_dir):
    files = []
    dir_prefix = working_dir + 'resources/battenberg/'
    # the running program
    files.append(dir_prefix+'beagle.22Jul22.46e.jar')
    # the folder containing all the reference files
    dir_b37 = dir_prefix + 'beagle/b37/'
    # the plink files
    chr_list = list(range(1,23))
    chr_list.extend(('X', 'X_par1', 'X_par2'))
    for chr in chr_list:
        chr = str(chr)
        # plink map files
        plink_map_file = dir_b37 + 'plink.chr' + chr + '.GRCh37.map'
        files.append(plink_map_file)
    # the bref3 files
    chr_list = chr_list[0:-2]
    for i in chr_list:
        i = str(i)
        bref_file = dir_b37 + 'chr' + i + '.1kg.phase3.v5a.b37.bref3'
        files.append(bref_file)
    return files

## 4. generate the output files after the preprocessing step of panConusig
def get_panConusig_preprocess(results, sampleID):
    files = []
    chr_list = list(range(1,22))
    for chr in chr_list:
        chr = str(chr)
        # essential output files from Battenberg
        alleleFreq = results + sampleID + '/06_panConusig/' + sampleID + '_alleleFrequencies_chr' + chr + '.txt'
        beagle_phase = results + sampleID + '/06_panConusig/' + sampleID + '_beagle5_output_chr' + chr + '.txt.vcf.gz'
        files.extend((alleleFreq, beagle_phase))
    # output from ASCAT.sc
    as_cna_raw = results + sampleID + '/06_panConusig/ASCAT_out/as_cna_profile_' + sampleID + '.sorted.dedup.bam_bam1.txt'
    as_cna_out = results + sampleID + '/06_panConusig/ASCAT_out/as_cna_profile_' + sampleID + '.txt'
    files.extend((as_cna_raw, as_cna_out))
    return files

## 5. generate the input file list for panConusig
def get_cna_profile(sample_df, results):
    files = []
    sample_list = sample_df.Sample
    for sampleID in sample_list:
        cna_profile = results + sampleID + '/06_panConusig/ASCAT_out/' + sampleID + '_as_cna_profile.tsv'
        files.append(cna_profile)
    return files












