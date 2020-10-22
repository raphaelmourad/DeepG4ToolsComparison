#Launch DeepG4 algorithm and compute score for each sequences (using package)
# author : Vincent ROCHER
# date : 05/10/20
if (!require(DeepG4))
{
  devtools::install_local("scripts/imports/DeepG4/",upgrade="never")
  if(!require(DeepG4)) stop("Package not found")
}
library(tidyverse)
library(DeepG4)
library(Biostrings)
args = commandArgs(trailingOnly=TRUE)

input_fas <- args[1]
output_file <- args[2]
custom_model <- args[3]
G4pos.seq=readDNAStringSet(input_fas)

resFreq <- Biostrings::letterFrequency(G4pos.seq, "N", as.prob = T)
testNFreq <- as.vector(resFreq > 0.1)
G4pos.seq <- G4pos.seq[!testNFreq,]


predictions <- DeepG4(G4pos.seq,model=custom_model)

tibble(names = names(G4pos.seq),score = predictions[,1]) %>% write_tsv(output_file)
