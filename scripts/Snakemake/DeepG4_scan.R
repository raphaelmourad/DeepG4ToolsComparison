#Launch DeepG4 algorithm and compute score for each sequences (using package)
# author : Vincent ROCHER
# date : 16/10/20
if (!require(DeepG4))
{
  devtools::install_local("scripts/imports/DeepG4/",upgrade="never")
  if(!require(DeepG4)) stop("Package not found")
}
library(tidyverse)
library(DeepG4)
library(Biostrings)
library(parallel)
args = commandArgs(trailingOnly=TRUE)

input_fas <- args[1]
output_file <- args[2]
custom_model <- args[3]
threads <- args[4]
G4pos.seq=readDNAStringSet(input_fas)


width.seq <- G4pos.seq %>% width %>% min
G4pos.seq.splitted <- mclapply(G4pos.seq,function(x){
  ExtractSubSequence(x = x, 
                     k = floor(width.seq/20))
},mc.cores=threads) %>% bind_rows(.id="names")

resFreq <- Biostrings::letterFrequency(DNAStringSet(G4pos.seq.splitted$seq), "N", as.prob = T)
testNFreq <- as.vector(resFreq > 0.1)
G4pos.seq.splitted <- G4pos.seq.splitted[!testNFreq,]

predictions <- DeepG4(DNAStringSet(G4pos.seq.splitted$seq),model=custom_model)

G4pos.seq.splitted <- G4pos.seq.splitted %>% 
  mutate(score = predictions[,1])  %>%
  group_by(names) %>% summarise(score = max(score)) %>% 
  write_tsv(output_file)
