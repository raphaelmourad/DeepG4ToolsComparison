# Use tensorflow model from penguinn to retrain with active G4
# Author : Vincent ROCHER
# Date : 14/10/2020
# From : scripts/generate_model_retrained/penquinn_retrained.R
library(keras)
library(reticulate)
library(tidyverse)
library(Biostrings)
source("scripts/functions/data_handling.R")
source_python("scripts/imports/penguinn/sequence_to_ohe.py")

args = commandArgs(trailingOnly=TRUE)


input_fas <- args[1]
ouput_file <- args[2]
model_path <- args[3]


#Load dataset
x_test <- as.character(readDNAStringSet(input_fas))
# Conversion into penguinn OneHot using sequence_to_ohe.py
x_test_array <- array(0,dim = c(length(x_test),200,4))
for(i in 1:length(x_test)){
  x_test_array[i,,] <- x_test[i]  %>% sequence_to_ohe
}

# load model
model <- load_model_hdf5(model_path)
# Prediction
res <- model %>% predict(x_test_array)
tibble("Sequence ID" = names(x_test),"Score"=res[,1]) %>% 
  write_tsv(ouput_file)
