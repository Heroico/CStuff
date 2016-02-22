#! /usr/bin/env Rscript
library(argparse)
source("Utilities.R")

parser <- ArgumentParser(description='Compare gene expressions')
parser$add_argument('--file',
                    help='Path of zscore',
                    default='results/YFS_GIANT.txt')


parser$add_argument('--result_prefix',
                    help='comparison prefix',
                    default='results/YFS_GIANT')

arguments <- parser$parse_args(commandArgs(TRUE))

input = arguments$file
output = paste(arguments$result_prefix, "qqunif.png", sep="", collapse="")

zscore_qqunif(input, output, " ")