#Vincent ROCHER (07/10/20)
#required packages
require(Biostrings)
require(itertools)
require(tidyverse)

args = commandArgs(trailingOnly=TRUE)


input_fas <- args[1]
input_table <- args[2]
output_table <- args[3]
type_calc <- args[4]

G4pos.seq=readDNAStringSet(input_fas)
G4pos.seq <- tibble(names = names(G4pos.seq))


G4pos.qparse <- read_tsv(input_table,col_names = F) %>%
  dplyr::select(X1,X3) %>% 
  group_by(X1) 


if(type_calc == "mean"){
  G4pos.qparse <- G4pos.qparse %>% summarise(score = mean(X3,na.rm = T)) 
}else if(type_calc == "max"){
  G4pos.qparse <- G4pos.qparse %>% summarise(score = max(X3,na.rm = T)) 
}else{
  G4pos.qparse <- G4pos.qparse %>% summarise(score = sum(X3,na.rm = T))
}


G4pos.qparse <- G4pos.seq %>% left_join(G4pos.qparse,by=c("names"="X1")) %>% replace_na(list(score = 0))


G4pos.qparse %>%  write_tsv(output_table)


