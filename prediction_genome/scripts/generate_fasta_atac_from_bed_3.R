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

getScoreBW <- function (WIG, BED,forScan=F)
{
  res <- do.call("rbind",lapply(split(BED, droplevels(GenomeInfoDb::seqnames(BED))), function(zz) {
    cov <- WIG[[unique(as.character(GenomeInfoDb::seqnames(zz)))]]
    score <- IRanges::Views(cov, start = BiocGenerics::start(zz), end = BiocGenerics::end(zz))
    return(as.matrix(score))
  }))
  if(forScan){
    return(res)
  }else{
    return(rowMeans(res))
  }
}


#read bed files
G4pos.bed <- snakemake@input[["genome"]] %>% read_bed()
my_seuil_bg <- snakemake@params[["seuil"]]
my_window_bg <- as.numeric(snakemake@params[["window"]])
selection <- IRanges::reduce(G4pos.bed)
#generate open chromatin dataset (ATAC or DNAse-seq)

ATAC_dataset <- snakemake@input[["atac_data"]]

G4pos.bed <- BiocGenerics::sort(GenomeInfoDb::sortSeqlevels(G4pos.bed))
G4all.atac <- ATAC_dataset %>% map(function(x){
  xbw <- x %>% import.bw(selection = BigWigSelection(selection),as="RleList")
  res <- xbw %>% getScoreBW(G4pos.bed)
  res[is.na(res)] <- 0
  return(res)
})%>% purrr::reduce(`+`) / length(ATAC_dataset)
message("step - get ATAC signal on bins: OK")
G4pos.bed.bg <- G4pos.bed %>% anchor_center() %>% mutate(width=my_window_bg)%>%
  as_tibble() %>%
  left_join(enframe(seqlens),by = c("seqnames"="name")) %>%
  mutate(end = ifelse(end>value,value,end)) %>% dplyr::select(-width) %>% as_granges()
G4all.atac.bg <- ATAC_dataset %>% map(function(x){
  xbw <- x %>% import.bw(selection = BigWigSelection(selection),as="RleList")
  res <- xbw %>% getScoreBW(G4pos.bed.bg)
  res[is.na(res)] <- 0
  return(res)
})%>% purrr::reduce(`+`) / length(ATAC_dataset)
message("step - get ATAC signal on background bins: OK")
my_test <- (G4all.atac/G4all.atac.bg)<my_seuil_bg
G4all.atac[my_test] <- 0
message("step - Replace values based on treshold: OK")
# saveRDS(G4all.atac,snakemake@output[["atac_merged"]])
tibble(atac_seq=G4all.atac) %>% write_tsv(snakemake@output[["atac_merged"]],col_names = F)
