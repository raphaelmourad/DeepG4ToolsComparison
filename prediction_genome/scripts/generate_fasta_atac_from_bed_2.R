# Vincent ROCHER
# created 05/03/2021
# CBI LBCMCP Legube Team, Toulouse
# create X fasta/atac file given a positive bed and control bed

#required packages
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


getScoreBW <- function (one.w, x,meanVal = T) 
{
  require(magrittr)
  res <- lapply(split(x, droplevels(seqnames(x))), function(zz) {
    message(unique(as.character(seqnames(zz))))
    cov <- one.w[[unique(as.character(seqnames(zz)))]]
    score <- IRanges::Views(cov, start = start(zz), end = end(zz))
    score %>% 
      as.matrix()
  }) %>% do.call("rbind",.)
  if(meanVal){
    return(rowMeans(res))
  }else{
    return(res)
  }
}

#read bed files
G4pos.bed <- snakemake@input[["genome"]] %>% read_bed()

selection <- IRanges::reduce(G4pos.bed)
#generate open chromatin dataset (ATAC or DNAse-seq)

ATAC_dataset <- snakemake@input[["atac_data"]]


G4all.atac <- ATAC_dataset %>% map(function(x){
  xbw <- x %>% import.bw(selection = BigWigSelection(selection),as="RleList")
  res <- xbw %>% getScoreBW(G4pos.bed)
  res[is.na(res)] <- 0
  return(res)
})%>% purrr::reduce(`+`) / length(ATAC_dataset)

saveRDS(G4all.atac,snakemake@output[["atac_merged"]])