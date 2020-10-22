# Matthieu Genais
# created on friday, 24/01/2020
# CBI LBCMCP Legube Team, Toulouse
# This code is made launch Quadron Snakemake rule
# Updated by Vincent ROCHER on 22/09/2020
#required packages
require(Biostrings)
require(itertools)
require(tidyverse)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)
#getting option for Quadron algorithm
type <- args[1]
Quadron_lib <- args[2]
source_lib <- args[3]
threads <- as.integer(args[4])
input_file <- args[5] 
output_file <- args[6]

source(source_lib)
registerDoParallel(cores=args[4])
#loading Fasta file
G4pos.seq=readDNAStringSet(input_file)
#G4pos.seq=readDNAStringSet("pipeline/data/samples/testElissar/Ct_test_Quadron.Fa")
#G4pos.seq=readDNAStringSet("/media/mourad/disk750Gb/G4pred_Matthieu/pipeline/results/201/Peaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b/fasta/merged/Peaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b_Ctrl_gkmSVM_merged.Fa")


#if the option is set to score
if(type=="score"){
  
  #then get the score
    
  Quadron_seq <- Quadron2(G4pos.seq,score = T,nCPU = threads)

}else{

  #get all features from Quadron output
  Quadron_seq <- Quadron2(G4pos.seq,score = F,nCPU = threads)
  colnames(Quadron_seq)=c("names",paste0("Quadron_",colnames(Quadron_seq)[2:length(colnames(Quadron_seq))]))


}

#export a TSV file for Snakemake
Quadron_seq %>%  write_tsv(output_file)