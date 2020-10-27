library(JASPAR2020)
library(PWMEnrich)
library(tidyverse)
library(TFBSTools)
library(ranger)
#From /media/ElissarDisk/Elissar/Projects/DeepG4/G4pred_Matthieu/scriptsR/Deep_learning_pipeline.R

rootMotifs = readMotifs("matrix-clustering_2020-04-30.103551_3iDIBq/matrix-clustering_cluster_root_motifs.tf",remove.acc = T)
namecluster <- names(rootMotifs)
rootMotifs <- lapply(1:length(rootMotifs),function(i){
    PFMatrix(ID=names(rootMotifs[i]),name =names(rootMotifs[i]) ,matrixClass = "PCM",profileMatrix = rootMotifs[[i]])
}) %>% c() %>% setNames(namecluster)

rootMotifs=do.call(TFBSTools::PFMatrixList,rootMotifs)

rootMotifsPFM = TFBSTools::toPWM(rootMotifs,type= 'prob')

# matrixList=TFBSTools::Matrix(rootMotifsPFM)
rdsfiles <- list.files("rds",pattern = "^Peaks_(BG4_G4seq|qG4-ChIP-seq-of-breast-cancer-PDTX)_.*_Sequence_train_test.rds",full.names = T)
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
    fileRF <- str_c("rds/RF_interpret_DeepG4_features_",basename(BG4dataset))
    if(file.exists(fileRF)){
        RF <- readRDS(fileRF)
    }else{
        RF=ranger(dependent.variable.name="y_train",data=dataRF,importance="impurity")
        saveRDS(RF,fileRF)
    }
    
    
    
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
    nbmotifs <- (colSums(mat_test)[match(names(impAll),colnames(mat_test))]/sum(mat_test))*100
    tibble(name = names(impAll),
              importanceVariable=impAll,
              abundanceKernel=nbmotifs,
           aucRF = as.numeric(aucRF))
},mc.cores=length(rdsfiles))

res <- res %>% setNames(str_extract(basename(rdsfiles),"BG4_G4seq_.+_GSE[0-9]+|qG4-ChIP-seq-of-breast-cancer-PDTX")) %>% bind_rows(.id="Peaks")


best_25_for_each <- res %>%
    group_by(Peaks) %>% arrange(desc(importanceVariable)) %>% dplyr::slice(1:25) %>% pull(name)

subres <- res %>%
    filter(name %in% best_25_for_each) %>%
    mutate(Peaks_AUC = str_c(Peaks," (",round(aucRF,3),")")) %>% 
    mutate(name = as.factor(name))
subres %>% 
    ggplot() +
    geom_point(aes(x=fct_reorder(name,importanceVariable,max),y=importanceVariable,size=abundanceKernel,col=name)) +
    geom_segment(aes(x=fct_reorder(name,importanceVariable,mean),xend=fct_reorder(name,importanceVariable,mean),y=0,yend=importanceVariable,col=name)) +
    coord_flip() +
    facet_wrap(~Peaks_AUC,scales="free_x",nrow=1) +
    theme_bw() +
    theme(legend.position = "none") -> p


library(ggseqlogo)
library(cowplot)
sequences <- readDNAStringSet(system.file("extdata", "test_G4_data.fa", package = "DeepG4"))

orderID <- subres %>% group_by(name) %>% summarise(cc=max(importanceVariable)) %>% arrange(desc(cc)) %>% pull(name)
p.pcm <- lapply(rootMotifs[orderID],function(x){as.matrix(x)})
p.pcm <- ggseqlogo(p.pcm,ncol=1) + theme_void(base_size=4) 
pf <- align_plots(p.pcm,p,align = "hv",axis="bt")
pdf("imp_variables_clusters_RF.pdf",height=24,width=18)
plot_grid(pf[[1]],pf[[2]],ncol=2,rel_widths = c(0.15,0.85))
dev.off()
#idxDF=data.frame(name=KernelWithoutMotif,idx=idxKernelWithoutMotif)

#df = left_join(df,idxDF,by=c("name"="name"))





plot = df %>% ggplot(aes(x=abundanceKernel,y=importanceVariable,label=name))+geom_point()+
    scale_x_continuous(trans = "log10")




