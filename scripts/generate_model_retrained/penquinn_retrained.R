# Use tensorflow model from penguinn to retrain with active G4
# Author : Vincent ROCHER
# Date : 23/09/2020
# From : /Documents/Vincent/DeepG4/penguinn/Retrain_penguinn.R
library(keras)
library(reticulate)
library(tidyverse)
library(Biostrings)
source("scripts/functions/data_handling.R")
source_python("scripts/imports/penguinn/sequence_to_ohe.py")

#Load dataset
BG4_data <- readRDS("rds/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_train_test.rds")
c(x_train, y_train) %<-% BG4_data$train
c(x_test, y_test) %<-% BG4_data$test
# Reshape sequence to 200 bp according to the model input (see penguinn paper)
x_train <- subseq(DNAStringSet(x_train),start=1,width=200) %>% as.character()
x_test <- subseq(DNAStringSet(x_test),start=1,width=200) %>% as.character()
# Conversion into penguinn OneHot using sequence_to_ohe.py
x_train_array <- array(0,dim = c(length(x_train),200,4))
for(i in 1:length(x_train)){
  x_train_array[i,,] <- x_train[i]  %>% sequence_to_ohe
}

#Load penguinn keras model architecture and retrained it using same hyperparameters (see penguinn paper)
my_config <- keras::get_config(load_model_hdf5("scripts/imports/penguinn/Models/model_1_1.h5"))
model <- keras::from_config(my_config)
model %>% compile(
  optimizer = optimizer_adam(lr=0.001, beta_1=0.9, beta_2=0.99),
  loss = "binary_crossentropy",
  metrics = c('accuracy')
)
#train the model
history <- model %>% fit(
  x_train_array,
  y_train,
  epochs = 15,
  batch_size = 32,
  validation_split = 0.2,
  verbose=1
)

model %>% save_model_hdf5("rds/penguinn_retrained_best_model_1_1.h5")
#Model evaluation

model <- load_model_hdf5("rds/penguinn_retrained_best_model_1_1.h5")
x_test_array <- array(0,dim = c(length(x_test),200,4))
for(i in 1:length(x_test)){
  x_test_array[i,,] <- x_test[i]  %>% sequence_to_ohe
}

c(pred_prob,pred_class) %<-% get_proba_classes(model,x_test_array,type = "Model")
c(table.1,p.1,p.2) %<-% evaluate_model(y_test,pred_prob,pred_class)
table.1
print(p.1)
print(p.2)
require(pROC)
rocG4detector <- roc(as.factor(y_test),pred_prob[,1],ci=T)
aucG4detector <- pROC::auc(rocG4detector)
# plot the roc and save it as pdf 
pdf(file=(paste0("AUCROC_penguinn_fine_tuning_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM.pdf")),width=3,height=3)
plot(rocG4detector,main = paste0("AUC = ",round(aucG4detector, digits = 3)))
dev.off()


