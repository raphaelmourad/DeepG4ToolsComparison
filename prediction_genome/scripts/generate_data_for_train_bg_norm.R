#Author: Vincent ROCHER
#Date: 19/04/2021

require(Biostrings)
if (!require(plyranges))
{
  BiocManager::install("plyranges")
}
require(plyranges)
if (!require(BSgenome.Hsapiens.UCSC.hg19))
{
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}
library(BSgenome.Hsapiens.UCSC.hg19)
require(tidyverse)
require(rsample)
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
my_seuil_bg <- snakemake@params[["seuil"]]
my_window_bg <- as.numeric(snakemake@params[["window"]])
dirG4 <- "../bed/"
oneG4 <- "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM.bed"
ATAC_dataset <- snakemake@input[["atac"]]

neg <- oneG4
pos <- str_replace(neg,"_Ctrl_gkmSVM.bed",".bed")

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

xbw <- ATAC_dataset %>% import.bw(as="RleList")
G4all.atac  <- c(
  xbw %>% getScoreBW(G4pos.bed ),
  xbw %>% getScoreBW(G4neg.bed )
)
G4all.atac[is.na(G4all.atac)] <- 0
G4pos.bed.bg <- G4pos.bed %>% anchor_center() %>% mutate(width=my_window_bg)
G4neg.bed.bg <- G4neg.bed %>% anchor_center() %>% mutate(width=my_window_bg)

G4all.atac.bg  <- c(
  xbw %>% getScoreBW(G4pos.bed.bg ),
  xbw %>% getScoreBW(G4pos.bed.bg )
)
G4all.atac.bg[is.na(G4all.atac.bg)] <- 0

my_test <- (G4all.atac/G4all.atac.bg)<my_seuil_bg
G4all.atac[my_test] <- 0
# G4neg.seq <- str_sub(G4neg.seq,1,201)

labelAll=c(rep(1,length(G4pos.seq)),rep(0,length(G4neg.seq)))
#Split Index
set.seed(my.seed)
SplitIndex <- G4all.seq %>% names() %>% ManualSplitIndex(percTrain = percTrain)
SplitIndex <- cbind(SplitIndex,label = labelAll[SplitIndex$name])
train.data <- SplitIndex %>% filter(Type == "training")
test.data <- SplitIndex %>% filter(Type == "test")

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
saveRDS(dataOneHot,snakemake@output[["rds"]])
