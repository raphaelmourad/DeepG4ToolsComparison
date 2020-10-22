#Launch G4iMGrinder algorithm and compute score for each sequences
# author : Vincent ROCHER
# date : 05/10/20
# Notes : installed from devtools::install_local("scripts/imports/G4iMGrinder/")
if (!require(G4iMGrinder))
{
  devtools::install_local("scripts/imports/G4iMGrinder/",upgrade ="never")
  if(!require(G4iMGrinder)) stop("Package not found")
}
library(G4iMGrinder)
library(tidyverse)
library(Biostrings)
args = commandArgs(trailingOnly=TRUE)

input_fas <- snakemake@input[["input_fas"]]
ouput_file <- snakemake@output[[1]]
threads <- snakemake@threads

G4pos.seq=readDNAStringSet(input_fas)
Rs  <- mclapply(G4pos.seq,function(seq){
    res <- G4iMGrinder(Name = "G4iMGrinder", Sequence = as.character(seq),G4hunter=F,PQSfinder = F,NCores = 1,Verborrea=F,Method3=F)
    if(nrow(res$PQSM2a) == 0){
        return(data.frame("freq"=NA,"Runs"=NA,"IL"=NA,"Sequence"=NA,"Length"=NA,"G"=NA,"T"=NA,"A"=NA,"C"=NA))
    }
    return(res$PQSM2a)
},mc.cores=threads)

Rs$PQSM2a %>% write_tsv(ouput_file)