# Generate Quadron features for Quadron retrained on active G4 dataset
# Script based on Elissar work on /media/ElissarDisk/Elissar/Projects/G4pred/Scripts/Scripts_R/AUROC/AUCROC_XGBoost_Quadron_trained.R
# Made by Vincent ROCHER (22/09/20)
library(keras) # To attribute value on multiple variables in one line
.libPaths(c("lib/Rlib",.libPaths()))
library(xgboost)
library(tidyverse)
library(Biostrings)

#Load dataset (already splited)
BG4_data <- readRDS("rds/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_train_test.rds")
c(x_train, y_train) %<-% BG4_data$train
c(x_test, y_test) %<-% BG4_data$test


Quadron_lib <- "scripts/imports/Quadron"
#Generate Quadron features
source("scripts/Snakemake/Quadron_add_Matthieu.R")
#Quadron train
Quadron.train <- Quadron2(DNAStringSet(x_train),score = F,nCPU = 4)
Quadron.train <- Quadron.train %>% mutate(Y = y_train) %>% drop_na()
saveRDS(Quadron.train,"rds/Quadron_features_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_train.rds")
#Quadron test
Quadron.test <- Quadron2(DNAStringSet(x_test),score = F,nCPU = 4)
Quadron.test <- Quadron.test %>% mutate(Y = y_test) %>% drop_na()
saveRDS(Quadron.test,"rds/Quadron_features_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_test.rds")
