# Take as input the bed format of G4CatchAll output and generate a tsv
# Author : Vincent ROCHER
#date 05/10/20
# Input bed file has the following columns:
#   1. description of the fasta sequence (e.g. NC_00024.11 Y chromosome)
# 2. start of the match
# 3. end of the match
# 4. size of the match
# 5. strand of the match (e.g. +)
# 6. positive strand sequence of the match (e.g. CCCTTCCCTTTCCCTCCC)
# 7. matched G-quadruplex-forming sequence (e.g. GGGAGGGAAAGGGAAGGG)
# 8. score of the matched G-quadruplex-forming sequence based on selected scoring scheme
library(tidyverse)
library(Biostrings)
args = commandArgs(trailingOnly=TRUE)


input_table <- args[1]
input_fas <- args[2]
ouput_file <- args[3]
type_calc <- args[4]

G4pos.seq=readDNAStringSet(input_fas)
G4pos.seq <- tibble(names = names(G4pos.seq))
G4pos.G4CatchAll <- read_tsv(input_table,col_names = F) %>% dplyr::select(X1,X8) %>%
  mutate(X8 = abs(X8)) %>% 
  group_by(X1)
if(type_calc == "mean"){
  G4pos.G4CatchAll <- G4pos.G4CatchAll %>% summarise(score = mean(X8,na.rm = T)) %>% replace_na(list(score = 0))
}else if(type_calc == "max"){
  G4pos.G4CatchAll <- G4pos.G4CatchAll %>% summarise(score = max(X8,na.rm = T)) %>% replace_na(list(score = 0))
}else{
  G4pos.G4CatchAll <- G4pos.G4CatchAll %>% summarise(score = sum(X8,na.rm = T))%>% replace_na(list(score = 0))
}

G4pos.G4CatchAll <- G4pos.seq %>% left_join(G4pos.G4CatchAll,by=c("names"="X1")) %>% replace_na(list(score = 0))

G4pos.G4CatchAll %>% write_tsv(ouput_file)