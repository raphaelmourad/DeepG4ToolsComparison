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
output_file <- args[2]
custom_model <- args[3]



G4pos.seq=readDNAStringSet(input_fas)
resFreq <- Biostrings::letterFrequency(G4pos.seq, "N", as.prob = T)
testNFreq <- as.vector(resFreq > 0.1)
my.names <- names(G4pos.seq)[!testNFreq]
G4pos.seq <- G4pos.seq[!testNFreq,] %>% DNAToNumerical()

model <- keras::load_model_hdf5(custom_model)
predictions <- model %>% predict(G4pos.seq)


tibble(names = my.names,score = predictions[,1]) %>% write_tsv(output_file)