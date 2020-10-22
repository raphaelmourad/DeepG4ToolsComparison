# Generate G4Hunter scores for G4Hunter Retrained on active G4 dataset
# Script based on scripts/Snakemake/G4Hunter.R
# Made by Vincent ROCHER (05/10/20)
modG4huntrefForEach <- function(seqi,hl=seq(from = 1, to = 2, by = 0.1)){
  sres <- lapply(set_names(hl,hl),function(i) modG4huntref(seqi,hl=i)) %>% map(as_tibble) %>% bind_rows(.id="hli") %>% 
    right_join(tibble(hli = as.character(hl)),by="hli") 
  return(sres)
}

require(keras)
require(tidyverse)
require(Biostrings)
require(GenomicRanges)

source("scripts/imports/G4Hunter/G4HunterRshinyscript.R")


#Load dataset (already splited)
BG4_data <- readRDS("rds/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_train_test.rds")
c(x_train, y_train) %<-% BG4_data$train
c(x_test, y_test) %<-% BG4_data$test

#Associate label to seq names 
seq_label_train <- tibble(seq = names(x_train),Y=y_train)
seq_label_test <- tibble(seq = names(x_test),Y=y_test)

#G4hunter parameters loading
start_threshold = 1
end_threshold = 2
pas = 0.1
#G4hunter is then launch
res.train <- mclapply(DNAStringSet(x_train),modG4huntrefForEach,hl=seq(from = start_threshold, to = end_threshold, by = pas),mc.cores=30)  %>% bind_rows(.id="seq") 
res.train <- res.train %>% 
  mutate(score = abs(score)) %>%
  replace_na(list(score = 0)) %>% dplyr::select(seq,hli,score) %>% 
  group_by(seq,hli) 

res.train <- res.train %>% summarise(score = max(score,na.rm = T)) 


res.train <- res.train %>%  spread(key=hli,value = score) %>% ungroup() %>% left_join(seq_label_train,by="seq")
saveRDS(res.train,"rds/G4Hunter_scores_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_train.rds")
#G4hunter test
res.test <- mclapply(DNAStringSet(x_test),modG4huntrefForEach,hl=seq(from = start_threshold, to = end_threshold, by = pas),mc.cores=30)  %>% bind_rows(.id="seq") 
res.test <- res.test %>% 
  mutate(score = abs(score)) %>%
  replace_na(list(score = 0)) %>% dplyr::select(seq,hli,score) %>% 
  group_by(seq,hli) 

res.test <- res.test %>% summarise(score = max(score,na.rm = T)) 


res.test <- res.test %>%  spread(key=hli,value = score) %>% ungroup() %>% left_join(seq_label_test,by="seq")
saveRDS(res.test,"rds/G4Hunter_scores_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_test.rds")
