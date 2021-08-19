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
#Functions (copied from DeepG4 package)
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
NormBW <- function(x,binbed){
  ranges_bin <-  base::range(IRanges::subsetByOverlaps(x,binbed)$score,na.rm = TRUE, finite = TRUE)
  x$score <- scales::rescale(x$score,to = c(0,1),from = ranges_bin)
  x <- IRanges::coverage(x,weight = "score")
  return(x)
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


G4all.bed <- c(G4pos.bed,G4neg.bed)
G4all.bed$order <- 1:length(G4all.bed)
G4all.bed <- BiocGenerics::sort(GenomeInfoDb::sortSeqlevels(G4all.bed))
G4all.atac <- ATAC_dataset %>% map(function(x){
  xbw <- x %>% import.bw %>% NormBW(binbed)
  res <- xbw %>% getScoreBW(G4all.bed)
  res[is.na(res)] <- 0
  return(res)
})%>% purrr::reduce(`+`) / length(ATAC_dataset)



G4all.bed.bg <- G4all.bed %>% anchor_center() %>% mutate(width=my_window_bg)%>%
  as_tibble() %>%
  left_join(enframe(seqlens),by = c("seqnames"="name")) %>%
  mutate(end = ifelse(end>value,value,end)) %>% dplyr::select(-width) %>% as_granges()

G4all.atac.bg <- ATAC_dataset %>% map(function(x){
  xbw <- x %>% import.bw() %>% NormBW(binbed)
  res <- xbw %>% getScoreBW(G4all.bed.bg)
  res[is.na(res)] <- 0
  return(res)
})%>% purrr::reduce(`+`) / length(ATAC_dataset)

my_test <- (G4all.atac/G4all.atac.bg)<my_seuil_bg
G4all.atac[my_test] <- 0


tibble(atac_seq=G4all.atac) %>% write_tsv(snakemake@output[["atac_merged"]],col_names = F)
