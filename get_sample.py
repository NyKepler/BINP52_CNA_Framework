# -*- coding: utf-8 -*-
#!/usr/bin/python3

"""
Title: get_sample.py
Author: Guyuan TANG
Date: 2023/11/06

Description: generate the sample tsv file from a mix metadata xlsx file. The tsv file of each type of samples will be generated separately.

List of packages/libraries:
    sys, pandas, numpy, os, re, glob

List of functions:
    sample_list(in_table, in_seq_dir, out_dir)


Steps:
1. generate a list to install all the unique types of samples
2. each unique type will create a sample tsv file
3. match the type with sample metadata to derive their libraries, use the libraries to seach for their reads
4. collect the paths and filenames of the reads, and output to sample tsv files 


Usage: Python get_sample.py --input [input_xlsx] --data [input_sequence_directory] --output [output_directory]
For help: Python get_sample.py --help


"""

import sys
import os
import pandas as pd
import numpy as np
import re
import glob




# Define a function to extract unique types of samples and search for the corresponding reads 
def sample_list(in_table, in_seq_dir, out_dir):
    df = pd.read_excel(in_table, header=0)
    type_list = df['Type'].unique().tolist()

    for t in type_list:
        out_file = out_dir + t + '.sample.tsv'
        with open(out_file, 'w') as out_file:
            # write the header
            print("sample_name\tpatient\tfastq_1\tfastq_2", file = out_file)
            # iterate by rows in the df to derive the sample_name and patient
            df_t = df[df['Type'] == t]
            for line in df_t.itertuples():
                sample_name = getattr(line, 'Library')
                patient = getattr(line, 'Patient')
                
                # generate the search patterns
                R1 = in_seq_dir + sample_name + '_*_R1_*.fastq.gz'
                R2 = in_seq_dir + sample_name + '_*_R2_*.fastq.gz'
                # get the path of the corresponding fastq files
                R1 = glob.glob(R1)[0]
                R2 = glob.glob(R2)[0]

                # print the information
                print('{}\t{}\t{}\t{}'.format(sample_name,patient,R1,R2), file = out_file)





# Recognize and take in all the input information
if len(sys.argv) == 2 and sys.argv[1] == '--help':
    print("\nUsage should be: Python get_sample.py --input [input_xlsx] --data [input_sequence_directory] --output [output_directory]\n")

elif ("--input" in sys.argv) and ("--data" in sys.argv) and ("--output" in sys.argv):
    # take in the input xlsx table
    in_loc = sys.argv.index("--input") + 1
    in_table = sys.argv[in_loc]

    # take in the directory storing the sequences
    in_seq_loc = sys.argv.index("--data") + 1
    in_seq_dir = sys.argv[in_seq_loc]
    if in_seq_dir[-1] != '/':
        in_seq_dir = in_seq_dir + '/'

    # take in the directory to store the output sample tsv files
    out_loc = sys.argv.index("--output") + 1
    out_dir = sys.argv[out_loc]
    if out_dir[-1] != '/':
        out_dir = out_dir + '/'

# Main script
    try:
        # to check whether the input files and directories are in correct type
        if os.path.isfile(in_table) and os.path.isdir(in_seq_dir) and os.path.isdir(out_dir):
            # run the main function
            sample_list(in_table, in_seq_dir, out_dir)

        else:
            print("Please use --help to check the usage again!")
    

    except FileNotFoundError as not_found:
        print("The file {} was not found!".format(not_found.filename))
    

