# -*- coding: utf-8 -*-
#!/usr/bin/python3

"""
Title: All_Sample_Solutions.py
Author: Guyuan TANG
Date: 2023/11/23

Description: generate the final table containing solutions of all the samples.

List of packages/libraries:
    sys, pandas, numpy, os

List of functions:
    table_generate(sample_table, method, max_solutions), sample_solutions(sample_dir, sample_table, method, max_solutions)


Steps:
1. take in thethe manual created sample table (copy from the table containing all the sample metadata) as a dataframe;
2. based on the method and maximum number of solutions, extend the sample table;
3. according to the sample types and names, search for the solution files;
4. derive the solution details (ploidy, cellularity, distance) and write to the sample table;
5. export and cover the sample table.


Usage: python All_Sample_Solutions.py --sample-list [sample_list_xlsx] --sample-dir [sample solutions directory] --method [method to calculate solutions] --max-solutions [maximum number of solutions] --binsize [binsize in kb]
For help: python All_Sample_Solutions.py --help

"""

import sys
import os
import pandas as pd
import numpy as np


# Define a function to extend columns in the sample table
def table_generate(sample_table, method, max_solutions, binsize):
    sample_table['BinSize'] = str(binsize) + 'kb'
    for i in range(1, max_solutions+1):
        new_ploidy = method + '_ploidy_' + str(i)
        new_cellularity = method + '_cellularity_' + str(i)
        new_MAD = method + '_MAD_' + str(i)

        sample_table = pd.concat([sample_table,pd.DataFrame(columns=[new_ploidy,new_cellularity,new_MAD])], sort=False)

    return sample_table


# Define a function to extract the solutions for each sample from the solution files
def sample_solutions(sample_dir, sample_table, method, max_solutions):
    for index,row in sample_table.iterrows():
        sample_name = index
        sample_type = row['Type']
        solution_file = sample_dir + sample_type + '/05_absolute_CN/' + sample_name + '/' + sample_name + '.solution.csv'
        solution = pd.read_csv(solution_file, header=0)
        if solution.shape[0] > 0:
            for n in range(solution.shape[0]):
                # when reach the maximum number of solutions
                if n >= max_solutions:
                    break
                else:
                    col_ploidy = method + '_ploidy_' + str(n+1)
                    col_cellularity = method + '_cellularity_' + str(n+1)
                    col_MAD = method + '_MAD_' + str(n+1)
                    # fill in the data
                    sample_table.loc[index, col_ploidy] = solution.loc[n, 'ploidy']
                    sample_table.loc[index, col_cellularity] = solution.loc[n, 'cellularity']
                    sample_table.loc[index, col_MAD] = solution.loc[n, 'distance']

    return sample_table



# Main script
if len(sys.argv) == 2 and sys.argv[1] == '--help':
    print("\nUsage should be: python All_Sample_Solutions.py --sample-list [sample_list_xlsx] --sample-dir [sample solutions directory] --method [method to calculate solutions] --max-solutions [maximum number of solutions] --binsize [binsize in kb]\n")

elif ("--sample-list" in sys.argv) and ("--sample-dir" in sys.argv) and ("--method" in sys.argv) and ("--max-solutions" in sys.argv) and ("--binsize" in sys.argv):
    # take in the sample list
    sample_table_loc = sys.argv.index("--sample-list") + 1
    sample_table_file = sys.argv[sample_table_loc]

    # take in the directory storing the solutions results
    dir_loc = sys.argv.index("--sample-dir") + 1
    sample_dir = sys.argv[dir_loc]
    if sample_dir[-1] != '/':
        sample_dir = sample_dir + '/'

    # take in the method used to calculate solutions (rascal or ichorCNA)
    method_loc = sys.argv.index("--method") + 1
    method = sys.argv[method_loc]

    # take in the maximum number of solutions
    max_solutions_loc = sys.argv.index("--max-solutions") + 1
    max_solutions = sys.argv[max_solutions_loc]

    # take in the binsize used in the workflow
    binsize_loc = sys.argv.index("--binsize") + 1
    binsize = sys.argv[binsize_loc]

    try:
        # check the inputs are all in correct forms
        if os.path.isfile(sample_table_file) and os.path.isdir(sample_dir):
            # generate the dataframe from sample table, with sample names as the indices
            sample_table = pd.read_excel(sample_table_file, header=0, index_col='Library')
            sample_table = table_generate(sample_table,method,max_solutions, binsize)
            # fill in the solutions for each sample
            sample_table = sample_solutions(sample_dir,sample_table, method, max_solutions)
            # export the final table
            sample_table.to_excel(sample_table_file)
        else:
            print("Please use --help to check the usage again!")
    except FileNotFoundError as not_found:
        print("The file {} was not found!".format(not_found.filename))



