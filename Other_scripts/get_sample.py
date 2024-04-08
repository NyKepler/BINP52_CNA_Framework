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
1. generate a list to contain all the unique types of samples
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
    
    out_file = out_dir + 'All_samples.tsv'
    with open(out_file, 'w') as out_file:
        # write the header
        print("sample_name\tpatient\ttype\tfastq_1\tfastq_2", file = out_file)
        # iterate by rows in the df to derive the sample_name and patient
        for line in df.itertuples():
            sample_name = getattr(line, 'Library')
            patient = getattr(line, 'Patient')
            sample_type = getattr(line, 'Type')
                
            for in_dir in in_seq_dir:
                # generate the search patterns
                R1 = in_dir + sample_name + '_*R1_*.fastq.gz'
                R2 = in_dir + sample_name + '_*R2_*.fastq.gz'
                    
                if glob.glob(R1) and glob.glob(R2): # to identify if the sample is in the directory
                    # get the path of the corresponding fastq files
                    R1 = glob.glob(R1)[0]
                    R2 = glob.glob(R2)[0]

                    # print the information
                    print('{}\t{}\t{}\t{}\t{}'.format(sample_name,patient,sample_type,R1,R2), file = out_file)




# Recognize and take in all the input information
if len(sys.argv) == 2 and sys.argv[1] == '--help':
    print("\nUsage should be: Python get_sample.py --input [input_xlsx] --data [input_sequence_directory, separated by comma if there are multiple input directories] --output [output_directory]\n")

elif ("--input" in sys.argv) and ("--data" in sys.argv) and ("--output" in sys.argv):
    # take in the input xlsx table
    in_loc = sys.argv.index("--input") + 1
    in_table = sys.argv[in_loc]

    # take in the directory storing the sequences
    in_seq_loc = sys.argv.index("--data") + 1
    in_seq_dir = sys.argv[in_seq_loc].split(',')
    for a in range(0, len(in_seq_dir)):
        if in_seq_dir[a][-1] != '/':
            in_seq_dir[a] = in_seq_dir[a] + '/'

    # take in the directory to store the output sample tsv files
    out_loc = sys.argv.index("--output") + 1
    out_dir = sys.argv[out_loc]
    if out_dir[-1] != '/':
        out_dir = out_dir + '/'

# Main script
    try:
        # to check whether the input files and directories are in correct type
        is_seq_dir = True
        for a in in_seq_dir:
            if not os.path.isdir(a):
                is_seq_dir = False

        if os.path.isfile(in_table) and os.path.isdir(out_dir) and is_seq_dir:
            # run the main function
            sample_list(in_table, in_seq_dir, out_dir)

        else:
            print("Please use --help to check the usage again!")
    

    except FileNotFoundError as not_found:
        print("The file {} was not found!".format(not_found.filename))
    