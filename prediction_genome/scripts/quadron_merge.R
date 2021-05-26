# Vincent ROCHER (29/03/2021)
# Merge quadron results (both score and retrained)
# required packages
require(Biostrings)
require(tidyverse)

input_table <- snakemake@input[["table"]]
output_file <- snakemake@output[[1]]

#read fasta file
input_table %>% map(read_tsv) %>% bind_rows() %>% write_tsv(output_file)