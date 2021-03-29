library(keras)
library(tidyverse)
library(Biostrings)
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

args = commandArgs(trailingOnly=TRUE)
input_fas <- args[1]
input_rds <- args[2]
output_file <- args[3]
custom_model <- args[4]



G4pos.seq=readDNAStringSet(input_fas)
ATAC.res <- readRDS(input_rds)
resFreq <- Biostrings::letterFrequency(G4pos.seq, "N", as.prob = T)
testNFreq <- as.vector(resFreq > 0.1)
my.names <- names(G4pos.seq)[!testNFreq]
G4pos.seq <- G4pos.seq[!testNFreq,] %>% DNAToNumerical()
ATAC.res <- ATAC.res[!testNFreq]

model <- keras::load_model_hdf5(custom_model)
predictions <- model %>% predict(list(G4pos.seq,ATAC.res))


tibble(names = my.names,score = predictions[,1]) %>% write_tsv(output_file)