# This code take as input the G4Hunter.R output and get mean/sum absolute G4Hunter score for each combination
# Of sequence*hl paramater
# Vincent ROCHER on 22/09/2020

#required packages
library(tidyverse)
library(Biostrings)
args = commandArgs(trailingOnly=TRUE)

input_table <- args[1]
ouput_file <- args[2]
type_calc <- args[3]
score_selec <- as.integer(args[4])

G4pos.G4Hunter <- read_tsv(input_table) %>% 
  mutate(score = abs(score)) %>% 
  replace_na(list(score = 0)) %>% dplyr::select(seq,hli,score) %>% 
  group_by(seq,hli) 


if(type_calc == "mean"){
  G4pos.G4Hunter <- G4pos.G4Hunter %>% summarise(score = mean(score,na.rm = T)) 
}else if(type_calc == "max"){
  G4pos.G4Hunter <- G4pos.G4Hunter %>% summarise(score = max(score,na.rm = T)) 
}else{
  G4pos.G4Hunter <- G4pos.G4Hunter %>% summarise(score = sum(score,na.rm = T))
}

G4pos.G4Hunter %>% dplyr::rename(names = seq) %>% filter(hli == score_selec) %>% dplyr::select(names,score)  %>% write_tsv(ouput_file)