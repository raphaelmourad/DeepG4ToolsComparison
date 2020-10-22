# Use tensorflow model from G4detector to retrain with active G4
# Author : Vincent ROCHER
# Date : 14/10/2020
# From : /Documents/Vincent/DeepG4/G4detector/retrained_G4detector.R
library(keras)
library(reticulate)
library(tidyverse)
library(Biostrings)
library(yardstick)
source("scripts/functions/data_handling.R")
source_python("scripts/imports/G4detector/oneHot.py")

#Load dataset
BG4_data <- readRDS("rds/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_train_test.rds")
c(x_train, y_train) %<-% BG4_data$train
c(x_test, y_test) %<-% BG4_data$test


x_train.num <- array(0,dim = c(length(x_train),1188,1))
for(i in 1:length(x_train)){
  x_train.num[i,,] <- x_train[[i]]  %>%oneHot() %>% array(dim=c(1188,1))
}

my_config <- keras::get_config(load_model_hdf5("scripts/imports/G4detector/models/random/model_rand_K_PDS.h5"))
model <- keras::from_config(my_config)
model %>% compile(
  optimizer = optimizer_adam(lr=1e-4, beta_1=0.9, beta_2=0.99, epsilon=1e-8, decay=1e-5),
  loss = "binary_crossentropy",
  metrics = c('accuracy')
)
#train the model
history <- model %>% fit(
  x_train.num,
  y_train,
  epochs = 15,
  batch_size = 128,
  validation_split = 0.2,
  verbose=1
)
model %>% save_model_hdf5("rds/G4detector_retrained_model_rand_K_PDS.h5")
#Model evaluation
model <- load_model_hdf5("rds/G4detector_retrained_model_rand_K_PDS.h5")
x_test.num <- array(0,dim = c(length(x_test),1188,1))
for(i in 1:length(x_test)){
  x_test.num[i,,] <- x_test[[i]]  %>%oneHot() %>% array(dim=c(1188,1))
}
source("scripts/functions/prediction_functions.R")
c(pred_prob,pred_class) %<-% get_proba_classes(model,x_test.num,type = "Model")
c(table.1,p.1,p.2) %<-% evaluate_model(y_test,pred_prob,pred_class)
table.1
print(p.1)
print(p.2)
require(pROC)
rocG4detector <- roc(as.factor(y_test),pred_prob[,1],ci=T)
aucG4detector <- pROC::auc(rocG4detector)
# plot the roc and save it as pdf 
pdf(file=snakemake@output[[1]],width=3,height=3)
plot(rocG4detector,main = paste0("AUC = ",round(aucG4detector, digits = 3)))
dev.off()

