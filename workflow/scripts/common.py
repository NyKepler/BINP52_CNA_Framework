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






