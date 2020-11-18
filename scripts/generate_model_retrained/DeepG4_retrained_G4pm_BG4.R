library(Biostrings)
library(DeepG4)
library(rsample)
library(keras)
# Read positive and negative set of sequences 

#Load dataset (already splited)
# BG4_data <- readRDS("rds/Peaks_G4seqpm_BG4_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_train_test.rds")
BG4_data <- readRDS("rds/Peaks_qG4_G4seq_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_train_test.rds")

c(x_train, y_train) %<-% BG4_data$train
c(x_test, y_test) %<-% BG4_data$test


training <- DeepG4(x_train,y_train,retrain=TRUE,retrain.path = "rds/DeepG4_retrained_qG4_G4seq.h5")

# test <- DeepG4(x_test,y_test,model = "rds/DeepG4_retrained_qG4_G4seq.h5")
test <- DeepG4(x_test,y_test,model = "rds/DeepG4_retrained_G4seqpm_BG4_1_1.h5")
library(cowplot)
p_res_train <- cowplot::plot_grid(plotlist = test[2:3])
print(p_res_train)
rocRFall <- pROC::roc(as.factor( y_test),test[[1]][,1],ci=T)
aucRF <- pROC::auc(rocRFall)

