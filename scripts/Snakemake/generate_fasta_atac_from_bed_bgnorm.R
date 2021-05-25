# Vincent ROCHER
# created 28/04/2021
# CBI LBCMCP Legube Team, Toulouse
# create X fasta/atac file given a positive bed and control bed
# adapted from generate_fasta_atac_from_bed.R
# ADD tresholding by background
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
binbed <- snakemake@params[["random_regions"]] %>% read_bed()
#Functions
NormBW <- function(x,binbed){
  ranges_bin <- x %>% filter_by_overlaps(binbed) %>% .$score  %>% range(na.rm = TRUE, finite = TRUE)
  x$score <- scales::rescale(x$score,to = c(0,1),from = ranges_bin)
  seqlengths(x) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(x)]
  x <- coverage(x,weight = "score")
  return(x)
}
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
readBed <- . %>% read_tsv(col_names = F) %>% dplyr::select(1:3) %>% setNames(c("seqnames","start","end")) %>% as_granges()
G4pos.bed <- snakemake@input[["bed_pos"]] %>% read_bed()
if(unique(width(G4pos.bed)) != 201){
  G4pos.bed <- snakemake@input[["bed_pos"]] %>% readBed
}
G4neg.bed <- snakemake@input[["bed_ctrl"]] %>% read_bed()
if(unique(width(G4neg.bed)) != 201){
  G4neg.bed <- snakemake@input[["bed_ctrl"]] %>% readBed
}

G4pos.bed <- G4pos.bed %>% filter(seqnames != "chrY")
G4neg.bed <- G4neg.bed %>% filter(seqnames != "chrY")

#generate open chromatin dataset (ATAC or DNAse-seq)
my_seuil_bg <- snakemake@params[["seuil"]]
my_window_bg <- as.numeric(snakemake@params[["window"]])
ATAC_dataset <- snakemake@input[["atac_data"]]

G4all.atac <- ATAC_dataset %>% map(function(x){
  xbw <- x %>% import.bw %>% NormBW(binbed)
  res <- c(
    xbw %>% getScoreBW(G4pos.bed),
    xbw %>% getScoreBW(G4neg.bed)
  )
  res[is.na(res)] <- 0
  return(res)
})%>% purrr::reduce(`+`) / length(ATAC_dataset)



G4pos.bed.bg <- G4pos.bed %>% anchor_center() %>% mutate(width=my_window_bg)%>%
  as_tibble() %>%
  left_join(enframe(seqlens),by = c("seqnames"="name")) %>%
  mutate(end = ifelse(end>value,value,end)) %>% dplyr::select(-width) %>% as_granges() 

G4neg.bed.bg <- G4neg.bed %>% anchor_center() %>% mutate(width=my_window_bg)%>%
  as_tibble() %>%
  left_join(enframe(seqlens),by = c("seqnames"="name")) %>%
  mutate(end = ifelse(end>value,value,end)) %>% dplyr::select(-width) %>% as_granges()

G4all.atac.bg <- ATAC_dataset %>% map(function(x){
  xbw <- x %>% import.bw() %>% NormBW(binbed)
  res <- c(
    xbw %>% getScoreBW(G4pos.bed.bg),
    xbw %>% getScoreBW(G4neg.bed.bg)
  )
  res[is.na(res)] <- 0
  return(res)
})%>% purrr::reduce(`+`) / length(ATAC_dataset)

my_test <- (G4all.atac/G4all.atac.bg)<my_seuil_bg
G4all.atac[my_test] <- 0

tibble(atac_seq=G4all.atac) %>% write_tsv(snakemake@output[["atac_merged"]],col_names = F)