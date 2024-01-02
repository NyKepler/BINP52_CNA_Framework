# Title: common.py
# Author: Guyuan TANG
# Date: 2023/12/29

# Description: the script to store functions used for generating the input and output file lists. Because some file names will have different variables (bin size) which are difficult to only use expand() function in Snakemake.

import pandas as pd

## 1. generate the final outputs for the whole workflow (part 2)
def get_output_all(sample_df, results):
    files = []
    # extract the unique groups from the sample_df
    groups = sample_df['Group'].unique().tolist()
    # each group with each signature method
    for groupID in groups:
        sub_df = sample_df[sample_df['Group'].isin([groupID])]
        binsize = sub_df['Binsize'].unique().tolist()
        binsize = binsize[0] + 'kb'
        # cnsignatures
        CN_sig = results + 'signatures/' + groupID + '/CNsig/' + groupID + '_CNsig.' + binsize + '.matrix.rds'
        files.append(CN_sig)
        # Pancan signatures
        PanSig = results + 'signatures/' + groupID + '/PanSig/' + groupID + '_PanSig.' + binsize + '.matrix.rds'
        files.append(PanSig)
    return files

## 2. generate the file list for output of absolute copy number profiles
def get_output_absolute(sample_df, results):
    files = []
    for sample_name in sample_df.Sample:
        binsize = sample_df.loc[sample_name,'Binsize']
        binsize = binsize + 'kb'
        file_name = results + sample_name + '/05_absolute_CN/' + sample_name + '_' + binsize + '_seg.tsv'
        files.append(file_name)
    return files

## 3. generate the input file list for each group in the signature validation steps
def get_group_file(sample_df, results, group):
    files = []
    sub_df = sample_df[sample_df['Group'].isin([group])]
    binsize = sub_df['Binsize'].unique().tolist()
    binsize = binsize[0] + 'kb'
    for sample_name in sub_df.Sample:
        file_name = results + sample_name + '/05_absolute_CN/' + sample_name + '_' + binsize + '_seg.tsv'
        files.append(file_name)
    return files




