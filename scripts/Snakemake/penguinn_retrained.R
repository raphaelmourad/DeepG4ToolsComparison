# Load a sequence and predict G4 based on penguinn retrained model
# Author : Vincent ROCHER
# Date : 23/09/20
# For snakemake

library(keras)
library(reticulate)
library(tidyverse)
library(Biostrings)
source("scripts/functions/prediction_functions.R")

model <- load_model_hdf5("rds/penguinn_retrained_best_model_1_1.h5")

x_test <- readDNAStringSet(snakemake@input[["fas"]])
y_test <- as.integer(grepl("pos",names(x_test)))
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
pdf(file=snakemake@output[[1]],width=3,height=3)
plot(rocG4detector,main = paste0("AUC = ",round(aucG4detector, digits = 3)))
dev.off()

