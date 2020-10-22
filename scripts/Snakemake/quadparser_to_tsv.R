# Matthieu Genais
# created on friday, 29/01/2020
# CBI LBCMCP Legube Team, Toulouse
# This code is made to launch quad_parser Snakemake rule
#Updated by Vincent ROCHER (07/10/20)
#required packages
require(Biostrings)
require(itertools)
require(tidyverse)

args = commandArgs(trailingOnly=TRUE)


input_table <- args[1]
input_fas <- args[2]
ouput_file <- args[3]
type_calc <- args[4]

#read fasta file
fas = readDNAStringSet(input_fas)


#read quadparser output
bed = read.table(input_table)

#count number of predicted G4  for all fasta sequences
# number_G4peaks = as.data.frame(sapply(names(fas), function(x,y)sum(str_count(y,x)),bed$V1))
quadparser_score <- tibble(names = names(fas)) %>% left_join(bed,by=c("names"="V1")) %>%
  dplyr::select(names,V5) %>% 
  dplyr::rename(score="V5") %>% 
  replace_na(list(score = 0)) %>% 
  group_by(names)


if(type_calc == "mean"){
  quadparser_score <- quadparser_score %>% summarise(score = mean(score,na.rm = T)) 
}else if(type_calc == "max"){
  quadparser_score <- quadparser_score %>% summarise(score = max(score,na.rm = T)) 
}else{
  quadparser_score <- quadparser_score %>% summarise(score = sum(score,na.rm = T))
}

#export a TSV file for Snakemake
quadparser_score %>%  write_tsv(ouput_file)




