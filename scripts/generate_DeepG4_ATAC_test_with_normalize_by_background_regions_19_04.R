#Author: Vincent ROCHER
#Date: 19/04/2021

library(keras)
library(tidyverse)
library(Biostrings)
library(rtracklayer)
library(rsample)
library(plyranges)
library(scales)
library(BSgenome.Hsapiens.UCSC.hg19)
# run once (control bin generations)
getSeq2 <- function(y,x){
  getSeq(x,y)
}
binbed <- read_bed("../bed/random_region_for_scaling_min_max.bed")
NormBW <- function(x,binbed){
  ranges_bin <- x %>% filter_by_overlaps(binbed) %>% .$score  %>% range(na.rm = TRUE, finite = TRUE)
  x$score <- scales::rescale(x$score,to = c(0,1),from = ranges_bin)
  seqlengths(x) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(x)]
  x <- coverage(x,weight = "score")
  return(x)
}
ManualSplitIndex <- function(Index=NULL,percTrain = 0.8,sampling=T){
  Index <- Index %>% enframe()
  if(sampling){
    Index <- Index[sample(1:nrow(Index)),]
  }
  Index <- Index %>% initial_split(prop = percTrain)
  
  return(
    list("training" = training(Index),
         "test" = testing(Index)) %>% bind_rows(.id = "Type") 
  )
}

getScoreBW <- function (one.w, x,meanVal = T) 
{
  
  if(meanVal){
    lapply(split(x, droplevels(seqnames(x))), function(zz) {
      message(unique(as.character(seqnames(zz))))
      cov <- one.w[[unique(as.character(seqnames(zz)))]]
      score <- IRanges::Views(cov, start = start(zz), end = end(zz))
      score %>%mean()
    }) %>% do.call("c",.)
  }else{
    lapply(split(x, droplevels(seqnames(x))), function(zz) {
      message(unique(as.character(seqnames(zz))))
      cov <- one.w[[unique(as.character(seqnames(zz)))]]
      score <- IRanges::Views(cov, start = start(zz), end = end(zz))
      score %>% 
        as.matrix()
    }) %>% do.call("rbind",.)
  }
}
DNAToNumerical <- function(x,tabv = c("T"=4,"G"=3,"C"=2,"A"=1),lower.case=F,seq.size = 201){
  if(lower.case){
    names(tabv) <- tolower(tabv)
  }
  x <- Biostrings::as.matrix(x)
  listMat <- list()
  for(i in 1:length(tabv)){
    nuc_index <- tabv[[i]]
    nuc_value <- names(tabv[i])
    mat <- matrix(0,nrow(x),ncol(x))
    mat[x==nuc_value] <- 1
    if(ncol(x)<seq.size){
      mat <- cbind(mat,matrix(0,nrow(x),seq.size-ncol(x)))
    }
    listMat[[nuc_index]] <- mat
  }
  arrayout <- array(unlist(listMat), dim = c(nrow(listMat[[1]]), ncol(listMat[[1]]), length(listMat)))
  return(arrayout)
  
}
#Parameters
my.seed <- 42
percTrain <- 0.8
my_seuil_bg <- 2
dirG4 <- "/home/disleg/Documents/Vincent/DeepG4/scripts/ATACseq_18012020/PIPELINE/data/BED/"
myG4 <- list.files(dirG4,pattern="Ctrl_gkmSVM.bed")
ATAC_seq_data <- list(
  "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM"=c("ATAC_entinostat_mean")
  # "Peaks_BG4_G4seq_HaCaT_GSE99205_hg19_201b_Ctrl_gkmSVM"=c("ATAC_entinostat_mean"),
  # "Peaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b_Ctrl_gkmSVM"=c(
  #   "rhh113_HEKnp_ATAC_701_517_24022015_normalized",
  #   "rhh114_HEKnp_ATAC_702_502_24022015_normalized",
  #   "rhh115_HEKnp_ATAC_703_503_24022015_normalized",
  #   "rhh116_HEKnp_ATAC_701_517_27032015_normalized",
  #   "rhh117_HEKnp_ATAC_702_502_27032015_normalized",
  #   "rhh118_HEKnp_ATAC_704_504_27032015_normalized"
  # ),
  # "Peaks_BG4_G4seq_K562_GSE107690_hg19_201b_Ctrl_gkmSVM"=c("GSM4133303_YET96_ATAC_K652"),
  # "Peaks_G4P_G4seq_GSE133379_293T_hg19_201b_Ctrl_gkmSVM"=c("ENCFF716SFD"),
  # # "Peaks_G4P_G4seq_GSE133379_A549_hg19_201b_Ctrl_gkmSVM"=c("SRX069099"),
  # "Peaks_G4P_G4seq_GSE133379_A549_hg19_201b_Ctrl_gkmSVM"=c("ENCFF180FXV"),
  # "Peaks_G4P_G4seq_GSE133379_H1975_hg19_201b_Ctrl_gkmSVM"=c(
  #   "GSM4217852_WT-rep1-ATAC",
  #   "GSM4217853_WT-rep2-ATAC"
  #   
  # ),
  # # "Peaks_G4P_G4seq_GSE133379_HeLaS3_hg19_201b_Ctrl_gkmSVM"=c("SRX100899"),
  # "Peaks_G4P_G4seq_GSE133379_HeLaS3_hg19_201b_Ctrl_gkmSVM"=c("SRX2370816"),
  # "Peaks_qG4_G4seq_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM"=c(
  #   "ENCFF615FRD",
  #   "ENCFF922TLC"
  # )
)
for(oneG4 in myG4){
  message(oneG4)
  ATAC_dataset <- str_c("../bigwig/",ATAC_seq_data[[str_remove(oneG4,".bed")]],".bw")
  
  neg <- oneG4
  pos <- str_replace(neg,"_Ctrl_gkmSVM.bed",".bed")
  # #Process output dir
  GenName <- neg %>% str_remove(".bed")
  output_dir <- "/home/disleg/Documents/Vincent/DeepG4/scripts/ATACseq_18012020/PIPELINE/data/rds/"
  output_prefix <- str_c(GenName,percTrain,my.seed,"rescale_BW_by_bg_5kb_seuil",my_seuil_bg,sep="_")
  
  #Load data
  readBed <- . %>% read_tsv(col_names = F) %>% dplyr::select(1:3) %>% setNames(c("seqnames","start","end")) %>% as_granges()
  G4pos.bed <- str_c(dirG4,pos) %>% read_bed()
  if(unique(width(G4pos.bed)) != 201){
    G4pos.bed <- str_c(dirG4,pos) %>% readBed
  }
  G4neg.bed <- str_c(dirG4,neg) %>% read_bed()
  if(unique(width(G4neg.bed)) != 201){
    G4neg.bed <- str_c(dirG4,neg) %>% readBed
  }
  G4pos.seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,G4pos.bed)
  names(G4pos.seq) <- str_c(str_remove(pos,".bed"),1:length(G4pos.seq),sep="_")
  G4neg.seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,G4neg.bed)
  names(G4neg.seq) <- str_c(str_remove(neg,".bed"),1:length(G4neg.seq),sep="_")
  G4all.seq=c(G4pos.seq,G4neg.seq)
  #Get ATAC-seq score for bed
  nb_cores <- ifelse(length(unlist(ATAC_dataset))>5,5,length(ATAC_dataset))
  
  
  G4all.atac <- ATAC_dataset %>% mclapply(function(x){
    xbw <- x %>% import.bw %>% NormBW(binbed)
    res <- c(
      xbw %>% getScoreBW(G4pos.bed ),
      xbw %>% getScoreBW(G4neg.bed )
    )
    res[is.na(res)] <- 0
    return(res)
  },mc.cores = nb_cores) %>% purrr::reduce(`+`) / length(ATAC_dataset)
  G4pos.bed.bg <- G4pos.bed %>% anchor_center() %>% mutate(width=5000)
  G4neg.bed.bg <- G4neg.bed %>% anchor_center() %>% mutate(width=5000)
  G4all.atac.bg <- ATAC_dataset %>% mclapply(function(x){
      xbw <- x %>% import.bw %>% NormBW(binbed)
      res <- c(
        xbw %>% getScoreBW(G4pos.bed.bg),
        xbw %>% getScoreBW(G4neg.bed.bg)
      )
      res[is.na(res)] <- 0
      return(res)
    },mc.cores = nb_cores) %>% purrr::reduce(`+`) / length(ATAC_dataset)
  
  
  my_test <- (G4all.atac/G4all.atac.bg)<my_seuil_bg
  G4all.atac[my_test] <- 0
  # G4neg.seq <- str_sub(G4neg.seq,1,201)
  
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
      "x_atac"= G4all.atac[train.data$name],
      "y"= train.data$label
    ),
    "test"=list(
      "x"= G4all.seq[test.data$name],
      "x_atac"= G4all.atac[test.data$name],
      "y"= test.data$label
    ))
  saveRDS(dataFasta,str_c(str_c(output_dir,"/Fasta/"),output_prefix,"_Sequence_train_test.rds"))
  
  
  G4all.num <- G4all.seq %>% DNAToNumerical
  
  dataOneHot <- list(
    "train"=list(
      "x"= G4all.num[train.data$name,,],
      "x_atac"= G4all.atac[train.data$name],
      "y"= train.data$label
    ),
    "test"=list(
      "x"= G4all.num[test.data$name,,],
      "x_atac"= G4all.atac[test.data$name],
      "y"= test.data$label
    ))
  saveRDS(dataOneHot,str_c(str_c(str_c(output_dir,"/oneHot/"),output_prefix,"_OneHot_train_test.rds")))
  
}
