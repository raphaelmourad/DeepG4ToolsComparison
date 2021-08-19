#Author: Vincent ROCHER
#Date: 19/04/2021
#Upated at 02/07/21 to use function from DeepG4 package
library(keras)
library(tidyverse)
library(Biostrings)
library(rtracklayer)
library(rsample)
library(plyranges)
library(scales)
library(BSgenome.Hsapiens.UCSC.hg19)

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




#Parameters
my.seed <- 42
percTrain <- 0.8
my_seuil_bg <- 2
dirG4 <- "bed/"
# myG4 <- list.files(dirG4,pattern="Ctrl_gkmSVM.bed")
myG4 <- "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM.bed"
ATAC_seq_data <- list(
  "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM"=c("ATAC_entinostat_mean")
)
message(oneG4)
ATAC_dataset <- str_c("bigwig/",ATAC_seq_data[[str_remove(oneG4,".bed")]],".bw")

neg <- oneG4
pos <- str_replace(neg,"_Ctrl_gkmSVM.bed",".bed")
# #Process output dir
GenName <- neg %>% str_remove(".bed")
output_dir <- "rds/"
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


ATAC <- import.bw(ATAC_dataset)
labelAll <- c(rep(1,length(G4pos.bed)),rep(0,length(G4neg.bed)))
c(G4all.seq,G4all.atac) %<-% DeepG4::DeepG4InputFromBED(BED = G4all.bed,ATAC = ATAC,is.bw = TRUE,GENOME = BSgenome.Hsapiens.UCSC.hg19,use.bg = TRUE,windows_bg=5000,treshold_bg = 2)

names(G4all.seq) <- c(str_c(str_remove(pos,".bed"),1:length(G4pos.bed),sep="_"),str_c(str_remove(neg,".bed"),1:length(G4neg.bed),sep="_"))

#Split Index
set.seed(my.seed)
SplitIndex <- G4all.seq %>% names() %>% ManualSplitIndex(percTrain = percTrain)
SplitIndex <- cbind(SplitIndex,label = labelAll[SplitIndex$name])
SplitIndex %>% write_tsv(str_c(output_dir,output_prefix,"_SplitIndex.tsv.gz"))
train.data <- SplitIndex %>% filter(Type == "training")
test.data <- SplitIndex %>% filter(Type == "test")

G4all.num <- G4all.seq %>% DeepG4::DNAToNumerical(tabv = c("T"=4,"G"=3,"C"=2,"A"=1),lower.case=F,seq.size = 201)

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
saveRDS(dataOneHot,str_c(str_c(output_dir,output_prefix,"_OneHot_train_test.rds")))

