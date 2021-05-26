# Vincent ROCHER (29/03/2021)
# Convert penguinn output into tsv format seq score for all chromosome
# required packages
require(Biostrings)
require(tidyverse)

input_table <- snakemake@input[["table"]]
output_file <- snakemake@output[[1]]

#read fasta file
input_table %>% map(read_tsv) %>% bind_rows() %>% write_tsv(output_file)