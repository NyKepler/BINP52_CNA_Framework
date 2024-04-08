#!/usr/bin/env Rscript

# Title: select_tab.R
# Author: Guyuan TANG
# Date: 2024/1/19

# Description: this script will be used to pick out the signature definition table from the Pan-Cancer study.

library(optparse)
library(openxlsx)

option_list <- list(
    make_option(c('-i','--input'), dest = "input_file",
                help = "the Excel file needed to be operated."),
    
    make_option(c('-o', '--output'), type = 'character', dest = "output_file", help = "the output filename.")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- "--help"

opt <- parse_args(option_parser, args)

input_file <- opt$input_file
output_file <- opt$output_file

if (is.null(input_file)) stop("The excel file needs to be specified.")

# load the supplementary table
def_df <- read.xlsx(input_file, sheet = 'ST_17_Signature_Comp_Defs', colNames = TRUE, rowNames = FALSE)

# output the table
write.xlsx(def_df, file = output_file)