# This code take as input the G4Hunter.R output and get probability to have a G4 based on G4Hunter ranger implementation
# Vincent ROCHER on 06/10/2020

#required packages
library(tidyverse)
library(Biostrings)
if (!require(ranger))
{
  install.packages("https://cran.r-project.org/src/contrib/ranger_0.12.1.tar.gz",repos=NULL,type="source")
  if(!require(ranger)) stop("Package not found")
}
library(ranger)




args = commandArgs(trailingOnly=TRUE)

model_path <- args[1]
input_table <- args[2]
ouput_file <- args[3]
type_calc <- args[4]

#Model loading
RangerBG4 <- readRDS(model_path)
#Data loading and preparation
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
#Detect label
G4pos.G4Hunter <- G4pos.G4Hunter %>%  spread(key=hli,value = score) %>% ungroup() 

seqnames <- G4pos.G4Hunter$seq

#proper column names for ranger
G4pos.G4Hunter <- G4pos.G4Hunter %>% dplyr::select(-seq) %>% as.data.frame()
colnames(G4pos.G4Hunter) <- paste("G4HunterScore",colnames(G4pos.G4Hunter),sep="_")
#Prediction
Test.Pred <- predict(RangerBG4, G4pos.G4Hunter)

tibble(names = seqnames,score = Test.Pred$predictions[,2]) %>% write_tsv(ouput_file)
