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


G4pos.seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,G4pos.bed)
names(G4pos.seq) <- paste0(str_remove(snakemake@input[["bed_pos"]],".bed"),"_",1:length(G4pos.seq),"_pos")
G4neg.seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,G4neg.bed)
names(G4neg.seq) <- paste0(str_remove(snakemake@input[["bed_pos"]],".bed"),"_",1:length(G4neg.seq),"_neg")

#merging the to fasta files
fasta_merged=c(G4pos.seq,G4neg.seq)
#output a Fasta file into Snakemake

writeXStringSet(fasta_merged,snakemake@output[["fasta_merged"]])

