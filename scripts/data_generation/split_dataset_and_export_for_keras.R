#Author: Vincent ROCHER
#Date: 04/02/2020
##Load fasta/control sequences of BG4 and generate input for Keras with sampling (train/test)
##UPDATED for paper in 22/09/2020
library(keras)
library(tidyverse)
library(Biostrings)
library(rtracklayer)
library(rsample)
#Functions
source("scripts/functions/data_handling.R")
#Parameters
my.seed <- 42
percTrain <- 0.8
kmer_size <- 3
dirG4 <- "fasta/"

pos <- "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b.Fa"
neg <- "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM.Fa"
# pos <- "Peaks_BG4_G4seq_HaCaT_GSE99205_hg19_201b.Fa"
# neg <- "Peaks_BG4_G4seq_HaCaT_GSE99205_hg19_201b_Ctrl_gkmSVM.Fa"
# pos <- "Peaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b.Fa"
# neg <- "PPeaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b_Ctrl_gkmSVM.Fa"
# pos <- "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b.Fa"
# neg <- "Peaks_BG4_G4seq_K562_GSE107690_hg19_201b_Ctrl_gkmSVM.Fa"
# #Process output dir
GenName <- neg %>% str_remove(".Fa")
output_dir <- "rds/"
output_prefix <- str_c(GenName,percTrain,my.seed,sep="_")

#Load data
G4pos.seq=as.character(readDNAStringSet(str_c(dirG4,pos)))
#for gkmsvm, sequence are 202 length (WHYYYYYY ?)
G4neg.seq=as.character(readDNAStringSet(str_c(dirG4,neg)))
# G4neg.seq <- str_sub(G4neg.seq,1,201)
names(G4neg.seq) <- str_c(names(G4neg.seq),"ctrl",sep="_")
G4all.seq=c(G4pos.seq,G4neg.seq)
labelAll=c(rep(1,length(G4pos.seq)),rep(0,length(G4neg.seq)))
#Split Index
set.seed(my.seed)
SplitIndex <- G4all.seq %>% names() %>% ManualSplitIndex(percTrain = percTrain)
SplitIndex <- cbind(SplitIndex,label = labelAll[SplitIndex$name])
SplitIndex %>% write_tsv(str_c(output_dir,output_prefix,"_SplitIndex.tsv.gz"))
train.data <- SplitIndex %>% filter(Type == "training")
test.data <- SplitIndex %>% filter(Type == "test")
#Output Sequence train/test
dataFasta <-  list(
  "train"=list(
    "x"= G4all.seq[train.data$name],
    "y"= train.data$label
  ),
  "test"=list(
    "x"= G4all.seq[test.data$name],
    "y"= test.data$label
  ))
saveRDS(dataFasta,str_c(output_dir,output_prefix,"_Sequence_train_test.rds"))
#Also write validation to test AUC on every tool on it
writeXStringSet(DNAStringSet(dataFasta$test$x[dataFasta$test$y == 1]),str_c("fasta/TestSet_",output_prefix,".Fa"))
writeXStringSet(DNAStringSet(dataFasta$test$x[dataFasta$test$y == 0]),str_c("fasta/TestSet_",output_prefix,"_Ctrl_neg.Fa"))

#Output for oneHot

G4all.num <- G4all.seq %>% DNAToNumerical %>% convertOneHot(tabv = c(1,2,3,4,5))

dataOneHot <- list(
  "train"=list(
    "x"= G4all.num[train.data$name,,],
    "y"= train.data$label
  ),
  "test"=list(
    "x"= G4all.num[test.data$name,,],
    "y"= test.data$label
  ))
saveRDS(dataOneHot,str_c(output_dir,output_prefix,"_OneHot_train_test.rds"))
