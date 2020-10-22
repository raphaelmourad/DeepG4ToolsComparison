# Matthieu Genais
# created on friday, 24/01/2020
# CBI LBCMCP Legube Team, Toulouse
# This code is made launch G4hunter Snakemake rule
# Updated by Vincent ROCHER on 22/09/2020

#Functions

#Loop over each hl which are all possible values for hl parameter in modG4huntref function
#return a tibble with G4Hunter algorithm for each given values of hl
#fill with NA when G4Hunter return no score -> no G4 where found
modG4huntrefForEach <- function(seqi,hl=seq(from = 1, to = 2, by = 0.1)){
  sres <- lapply(set_names(hl,hl),function(i) modG4huntref(chr=seqi,hl=i)) %>% map(as_tibble) %>% bind_rows(.id="hli") %>% 
    right_join(tibble(hli = as.character(hl)),by="hli") 
  return(sres)
}

#required packages

require(tidyverse)
require(Biostrings)
require(GenomicRanges)

source(snakemake@params[["source"]])


#fasta loading
seq <- readDNAStringSet(snakemake@input[["fas"]]) 

#G4hunter parameters loading
x <-  snakemake@params[["start_threshold"]] %>% as.numeric()
y <-  snakemake@params[["end_threshold"]] %>% as.numeric()
z <-  snakemake@params[["pas"]] %>% as.numeric()
threads <- as.integer(snakemake@threads)
#G4hunter is then launch
res <- mclapply(seq,modG4huntrefForEach,hl=seq(from = x, to = y, by = z),mc.cores=threads)  %>% bind_rows(.id="seq") 
res %>% write_tsv(snakemake@output[["out"]])


