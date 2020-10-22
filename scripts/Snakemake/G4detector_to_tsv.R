# Vincent ROCHER (07/10/20)
# Convert G4 detector output into tsv format seq score
# required packages
require(Biostrings)
require(tidyverse)

args = commandArgs(trailingOnly=TRUE)


input_table <- args[1]
input_fas <- args[2]
ouput_file <- args[3]


#read fasta file
fas = readDNAStringSet(input_fas)

G4pos.G4detector <- read_delim(input_table,delim = ",")

tibble(names = names(fas),seq = as.character(fas)) %>% 
  left_join(G4pos.G4detector,by ="seq") %>% dplyr::select(names,score) %>% 
  replace_na(list(score = 0)) %>% 
  write_tsv(ouput_file)