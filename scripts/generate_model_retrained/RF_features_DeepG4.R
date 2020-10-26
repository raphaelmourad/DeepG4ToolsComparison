library(JASPAR2020)
library(PWMEnrich)
library(tidyverse)
library(TFBSTools)
library(ranger)
#From /media/ElissarDisk/Elissar/Projects/DeepG4/G4pred_Matthieu/scriptsR/Deep_learning_pipeline.R

rootMotifs = readMotifs("matrix-clustering_2020-04-30.103551_3iDIBq/matrix-clustering_cluster_root_motifs.tf",remove.acc = T)

rootMotifs <- lapply(1:length(rootMotifs),function(i){
    PFMatrix(ID=names(rootMotifs[i]),name =names(rootMotifs[i]) ,matrixClass = "PCM",profileMatrix = rootMotifs[[i]])
}) %>% c()

rootMotifs=do.call(TFBSTools::PFMatrixList,rootMotifs)

rootMotifsPFM = TFBSTools::toPWM(rootMotifs,type= 'prob')

# matrixList=TFBSTools::Matrix(rootMotifsPFM)
rdsfiles <- list.files("rds",pattern = "Peaks_BG4_G4seq_.*_Sequence_train_test.rds",full.names = T)
res <- mclapply(rdsfiles,function(BG4dataset){
    
    #Load dataset (already splited)
    BG4_data <- readRDS(BG4dataset)
    x_train <- BG4_data$train$x
    y_train <- BG4_data$train$y
    
    #matchMotifs function
    ix = motifmatchr::matchMotifs(rootMotifsPFM,DNAStringSet(x_train),out = "scores")
    
    #count number of motifs in each sequence
    mat = motifmatchr::motifCounts(ix)
    
    #change into dataframe
    dataKernel = as.data.frame(as.matrix(mat))
    colnames(dataKernel) = ID(rootMotifsPFM)
    
    dataKernel <- dataKernel %>% as("Matrix")
    
    
    # Random forests
    dataRF=cbind(y_train,dataKernel)
    RF=ranger(dependent.variable.name="y_train",data=dataRF,importance="permutation")
    saveRDS(RF,str_c("rds/RF_interpret_DeepG4_features_",BG4dataset))
    
    
    # Predict
    x_test <- BG4_data$test$x
    y_test <- BG4_data$test$y
    #count number of motifs in each sequence
    mat_test <- motifmatchr::motifCounts(motifmatchr::matchMotifs(rootMotifsPFM,DNAStringSet(x_test),out = "scores"))
    
    #change into dataframe
    mat_test = as.data.frame(as.matrix(mat_test))
    colnames(mat_test) = ID(rootMotifsPFM)
    
    mat_test <- mat_test %>% as("Matrix")
    
    ypred_prob= predict(RF,data=mat_test)$predictions
    rocRFall <- pROC::roc(as.factor( y_test),ypred_prob,ci=T)
    aucRF <- pROC::auc(rocRFall)
    
    impAll <- RF$variable.importance
    
    tibble(name = names(impAll),
              importanceVariable=impAll,
              abundanceKernel=colSums(dataKernel)[match(names(impAll),colnames(dataKernel))],aucRF = aucRF)
},mc.cores=length(rdsfiles))

#idxDF=data.frame(name=KernelWithoutMotif,idx=idxKernelWithoutMotif)

#df = left_join(df,idxDF,by=c("name"="name"))





plot = df %>% ggplot(aes(x=abundanceKernel,y=importanceVariable,label=name))+geom_point()+
    scale_x_continuous(trans = "log10")




