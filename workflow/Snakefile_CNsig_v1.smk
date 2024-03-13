# Title: Snakefile_CNsig.smk
# Author: Guyuan TANG
# Date: 2023/12/27 - 2024/1/2

import pandas as pd

from scripts.common import get_output_absolute, get_ref_battenberg, get_ref_beagle, get_panConusig_preprocess, get_cna_profile

# specify the configuration file
configfile: "config/config.yaml"

# specify the working directory
working_dir = config['workdir']
if working_dir[-1] != '/':
    working_dir = working_dir + '/'

# specify the results location (output directory)
results = config['outputdir']
if results[-1] != '/':
    results = results + '/'

# specify the sample information
sample_df = (pd.read_csv(config['samp_solutions'], 
    sep='\t',
    dtype={'Sample':str, 'Patient':str, 'Type':str, 'Group':str, 'Binsize':str, 'rds':str, 'Ploidy':float, 'Cellularity':float})
    .set_index('Sample', drop=False))

# specify the final output of the workflow
"""
The final output of this workflow should be the signature similarity matrix for each type of signature.
"""
rule all:
    input:
        CN_sig_SS = results + 'signatures/CN_sig/CN_sig.SSmatrix.rds',
        PanCan_sig_SS = results + 'signatures/PanCan_sig/PanCan_sig.SSmatrix.rds',
        panConusig_SS = results + 'signatures/panConusig/panConusig.SSmatrix.rds'


########## 1 Absolute copy number profiles ####################
"""
The final output for this step should be a segmented copy number profile (tsv) for each sample.
"""
rule absolute_seg_CN:
    input:
        files = get_output_absolute(sample_df, results)

# 1.1 prepare the rascal environment
rule rascal_env:
    output:
        "log/rascal_settle_info2.txt"
    conda: 'envs/rascal.yaml'
    script: 'scripts/rascal_env.R'

# 1.2 calculate the segmented absolute copy number profiles for each sample
rule rascal_absolute_CN:
    ### rascal will be used to calculate the absolute copy numbers and to apply segmentation
    input:
        env_set = 'log/rascal_settle_info2.txt',
        sample_table = config['samp_solutions'],
        rds = lambda wildcards: sample_df.loc[wildcards.sample, 'rds']
    output:
        results + '{sample}/05_absolute_CN/{sample}_{binsize}kb_seg.tsv'
    wildcard_constraints:
        binsize="\d+"
    params:
        sample = '{sample}',
        outdir = results + '{sample}/05_absolute_CN/'
    threads: 5
    conda: 'envs/rascal.yaml'
    script: 'scripts/fit_absoluteCN.R'


########## 2 CN Signatures ####################
"""
The final output for this step would be the matrix files (including a matrix txt file, a matrix object RDS, a simple heatmap for both sample-by-component matrix and sample-by-signature matrix) for each group. The matrices containing sample-by-signature information.
"""
# 2.1 clone the github repository used for signature validation
rule cn_sig_git:
    output:
        'workflow/scripts/cnsignatures/main_functions.R',
        'workflow/scripts/cnsignatures/helper_functions.R'
    shell: '''
    git clone https://bitbucket.org/britroc/cnsignatures.git 
    mv cnsignatures/ workflow/scripts/
    '''

# 2.2 check all the files are well prepared
rule cn_sig_check:
    input:
        link_up = rules.absolute_seg_CN.input,
        scripts = rules.cn_sig_git.output
    output:
        'log/cn_sig_settle.txt'
    shell: """
    echo 'Finished preparation for signature validation.' > {output}
    """

# 2.3 validate the signatures in our samples
rule CN_signature:
    input:
        # check whether the essential git repository have been downloaded
        main = 'workflow/scripts/cnsignatures/main_functions.R',
        helper = 'workflow/scripts/cnsignatures/helper_functions.R',
        # check all the absolute copy number profiles have successfully been generated
        check_data = rules.cn_sig_check.output,
        # the input segment files
        seg_files = rules.absolute_seg_CN.input
    output:
        CN_sig = results + 'signatures/CN_sig/CN_sig.SSmatrix.rds'
    params:
        sample_info = config['samp_solutions'],
        def_SC = 'workflow/scripts/cnsignatures/data/feat_sig_mat.rds',
        indir = results,
        outdir = results + 'signatures/CN_sig/'
    threads: 20
    conda: 'envs/CNsig.yaml'
    script: 'scripts/CN_sig.R'


########## 3 Pan-Cancer Signatures ####################
"""
The final output for this step should be the validated Pan-cancer signature matrix (sample-by-component) for each group. The actural outputs include a sample-by-component matrix txt, a full object RDS containing all information (such as activities and weights), two heatmaps (one for activities and one for sample-by-component).
"""
# 3.1 download the PanCan signature definiation table
rule PanCan_def_tab:
    # user should download the supplements manually
    # website: https://www-nature-com.ludwig.lub.lu.se/articles/s41586-022-04789-9#Sec17
    input:
        config['supple_tab']
    output:
        'resources/PanCan.xlsx'
    params:
        select_tab = 'workflow/scripts/select_tab.R',
        folder = config['supple_folder']
    conda: 'envs/PanCanSig.yaml'
    shell:'''
        mv '{input}' resources/Supplementary_Tables_15-22.xlsx
        Rscript {select_tab} -i resources/Supplementary_Tables_15-22.xlsx -o resources/PanCan.xlsx
        rm -r {folder}
    '''

# 3.2 set up the snakemake environment for running the validation
rule PanCan_sig_env:
    input:
        'resources/PanCan.xlsx'
    output:
        "log/PanCan_settle_info.txt"
    conda: 'envs/PanCanSig.yaml'
    script: 'scripts/PanCan_sig_env.R'

# 3.3 use the CINSignatureQuantification package to validate the signatures
rule PanCan_sig:
    input:
        env_set = "log/PanCan_settle_info.txt",
        seg_files = rules.absolute_seg_CN.input
    output:
        PanCan_sig_SS = results + 'signatures/PanCan_sig/PanCan_sig.SSmatrix.rds'
    params:
        sample_info = config['samp_solutions'],
        def_SC = 'resources/PanCan.xlsx',
        indir = results,
        outdir = results + 'signatures/PanCan_sig/'
    threads: 10
    conda: 'envs/PanCanSig.yaml'
    script: 'scripts/PanCan_sig.R'


########## 4 panConusig ####################
"""
The final output for this part should be the sample-by-component matrix for panConusig.
"""
# 4.1 prepare the environment (for Battenberg, ASCAT.sc and panConusig packages)
rule panConusig_env:
    output:
        "log/panConusig_settle_info.txt"
    conda: 'envs/panConusig.yaml'
    script: 'scripts/panConusig_env.R'

# 4.2 prepare the reference files for Battenberg and the tool beagle5
rule panConusig_ref_prep:
    input:
        env_set = "log/panConusig_settle_info.txt",
        impute_00 = working_dir + 'resources/battenberg/battenberg_impute_v3/impute_info00.txt'
    output:
        ref_battenberg = get_ref_battenberg(working_dir),
        ref_beagle = get_ref_beagle(working_dir)
    params:
        workdir = working_dir,
        impute_info = working_dir + 'resources/battenberg/impute_info.txt'
    shell: '''
    cat {input.impute_00} | sed 's#<path_to_impute_reference_files>#{params.workdir}#g' > {params.impute_info}
    '''

# 4.3 preprocessing steps (Battenberg and ASCAT.sc) before panConusig
rule panConusig_preprocessing:
    input: 
        ref_files = rules.panConusig_ref_prep.output,
        usr_battenberg = "workflow/scripts/usr_battenberg.R",
        bam_file = results + '{sample}/03_clean_up/{sample}.sorted.dedup.bam'
    output:
        as_cna_out = results + '{sample}/06_panConusig/ASCAT_out/{sample}_as_cna_profile.tsv'
    params:
        sampleID = '{sample}',
        impute_ref_dir = working_dir + 'resources/battenberg',
        beagle_ref_dir = working_dir + 'resources/battenberg/beagle',
        result_dir = results + '{sample}/06_panConusig/',
        ASCAT_outdir = results + '{sample}/06_panConusig/ASCAT_out/',
        binsize = lambda wildcards: sample_df.loc[wildcards.sample, 'Binsize']
    threads: config['battenberg_threads']
    conda: 'envs/panConusig.yaml'
    script: 'scripts/panConusig_preprocess.R'

# 4.4 generating the panConusig
rule panConusig_sig:
    input:
        cna_profiles = get_cna_profile(sample_df, results)
    output:
        panConusig_SS = results + 'signatures/panConusig/panConusig.SSmatrix.rds'
    params:
        def_SC = 'resources/panConusig_id.txt',
        outdir = results + 'signatures/panConusig/'
    conda: 'envs/panConusig.yaml'
    script: 'scripts/panConusig.R'




