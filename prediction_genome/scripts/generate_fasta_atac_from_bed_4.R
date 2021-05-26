# Vincent ROCHER
# created 19/04/2021
# CBI LBCMCP Legube Team, Toulouse
# create X fasta/atac file given a positive bed and control bed
#The idea here is to avoid background signal, so we set the accessibility signal to 0 if the ratio between peak/(peak_5kb) is below 2
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

seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
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

#read bed files
G4pos.bed <- snakemake@input[["genome"]] %>% read_bed()
my_seuil_bg <- 2
selection <- IRanges::reduce(G4pos.bed)
#generate open chromatin dataset (ATAC or DNAse-seq)

ATAC_dataset <- snakemake@input[["atac_data"]]
peaks_dataset <- snakemake@input[["atac_peaks"]]

G4all.atac <- ATAC_dataset %>% map(function(x){
  xbw <- x %>% import.bw(selection = BigWigSelection(selection),as="RleList")
  res <- xbw %>% getScoreBW(G4pos.bed)
  res[is.na(res)] <- 0
  return(res)
})%>% purrr::reduce(`+`) / length(ATAC_dataset)


#Peaks percentage coverage
peaks_coverage <- tryCatch(
  {
    read_bed(peaks_dataset)%>% coverage()
  },
  error=function(cond) {
    read_narrowpeaks(peaks_dataset)%>% coverage()
  }
)

G4all.peaks_coverage <- peaks_coverage %>% getScoreBW(G4pos.bed )

G4all.peaks_coverage[is.na(G4all.peaks_coverage)] <- 0

saveRDS(list(G4all.atac,G4all.peaks_coverage),snakemake@output[["atac_merged"]])
