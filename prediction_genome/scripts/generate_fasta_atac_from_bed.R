# Vincent ROCHER
# created 05/03/2021
# CBI LBCMCP Legube Team, Toulouse
# create X fasta/atac file given a positive bed and control bed

#required packages

require(tidyverse)

args = commandArgs(trailingOnly=TRUE)

input_tab <- args[1]
output_file <- args[2]


read_tsv(input_tab) %>% dplyr::select(-1:-3) %>% rowMeans() %>% saveRDS(output_file)