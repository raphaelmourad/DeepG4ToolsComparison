# Elissar Nassereddine 
# created on Tuesday, 14/11/2019
# CBI LBCMCP Legube Team, Toulouse
# This code load the original peaks for + and - reads of G4seq to merge them into one granges bed file as output 



setwd(dir = "/media/ElissarDisk/Elissar/Projects/DeepG4/Pipeline_Data/")

pacman::p_load(BSgenome.Hsapiens.UCSC.hg19,plyranges,GenomicRanges,readr,tidyr,tidyverse,dplyr)
genome <- BSgenome.Hsapiens.UCSC.hg19
win <- 100


#### to filtrate bed according to chromosome names 
Seqnames<- stringr::str_c("chr",c(1:22,"X"))

## list the originial files of plus and minus peaks of G4seq
PeaksPlusMinus <- list.files(paste0("Orig_Data/G4seq"),recursive=T,full.names=T,pattern = "intersect.bed.gz")

## read both files, concatenate them, set the column names and get the granges,finally reduce the granges(get rid of comon peaks) 
PeaksG4seq.GR= PeaksPlusMinus %>% purrr::map(read_tsv,col_names = F)%>% dplyr::bind_rows() %>% 
    setNames(c("seqnames","start","end")) %>% GRanges() %>% GenomicRanges::reduce()

##  get rid of unwanted chromosomes 
PeaksG4seq.GR.Filtered <- PeaksG4seq.GR[(seqnames(PeaksG4seq.GR) %in% Seqnames)]
as_tibble(PeaksG4seq.GR.Filtered)[,1:3] %>% write_tsv("Bed/Peaks_G4seq_hg19.bed", col_names = FALSE)


# Extract fasta from bed coordinates
Peaks.seq <- getSeq(genome,PeaksG4seq.GR.Filtered)
names(Peaks.seq)=paste0("G4seq_peak_",1:length(Peaks.seq))

# write fasta files
writeXStringSet(Peaks.seq,"Fasta/Peaks_G4seq_hg19.Fa")



  PeaksG4seq.GR.201 <-  PeaksG4seq.GR.Filtered %>% anchor_center() %>% mutate(width = win*2+1) 
  
  as_tibble(PeaksG4seq.GR.201)[,1:3] %>% write_tsv(paste0("Bed/Peaks_G4seq_hg19_",win*2+1,"b.bed"), col_names = FALSE)

  
  Peaks.seq.201 <- getSeq(genome,PeaksG4seq.GR.201)
  names(Peaks.seq.201)=paste0("G4seq_peak_",win*2+1,"b_",1:length(Peaks.seq.201))
  
  
  # write fasta files
  writeXStringSet(Peaks.seq.201,paste0("Fasta/Peaks_G4seq_hg19_",win*2+1,"b.Fa"))

  

