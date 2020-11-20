
# Elissar Nassereddine 
# created on Tuesday,26/10/2020
# CBI LBCMCP Legube Team, Toulouse
# This code load the original peaks for + and - reads of G4seq to merge them into 
# one granges bed and fasta files as output 



setwd(dir = "/media/ElissarDisk/Elissar/Projects/DeepG4/Pipeline_Data/")
pacman::p_load(BSgenome.Hsapiens.UCSC.hg19,plyranges,GenomicRanges,readr,tibble)
genome <- BSgenome.Hsapiens.UCSC.hg19

win = 100
#### to filtrate bed according to chromosome names 
Seqnames<- stringr::str_c("chr",c(1:22,"X"))


## list the originial files of plus and minus peaks of G4seq
PeaksPlusMinus <- list.files(paste0("Orig_Data/G4seq"),recursive=T,full.names=T,pattern = "intersect.bed.gz")

## read both files, concatenate them, set the column names and get the granges,finally reduce the granges(get rid of comon peaks) 
PeaksG4seq.GR= PeaksPlusMinus %>% purrr::map(read_tsv,col_names = F)%>% dplyr::bind_rows() %>% 
  purrr::set_names(c("seqnames","start","end")) %>% plyranges::as_granges()


##  get rid of unwanted chromosomes 
PeaksG4seq.GR.Filtered <- PeaksG4seq.GR[(seqnames(PeaksG4seq.GR) %in% Seqnames)]

## TO generate Peaks of 201 bp 
## use anchor center 
Peaks_GR <-  PeaksG4seq.GR.Filtered %>% anchor_center() %>% mutate(width = (win*2)+1) 
Peaks_201 <- Peaks_GR  %>% as_tibble() %>% dplyr::select(seqnames,start,end) %>%  write_tsv("Bed/Peaks_G4seqpm_hg19_201b.bed", col_names = FALSE)



### Generate fastas 
G4seq201 <- getSeq(genome,Peaks_GR) 
names(G4seq201) <- paste0("Peaks_G4seqpm_",1:length(G4seq201))
writeXStringSet(G4seq201,paste0("Fasta/Peaks_G4seqpm_hg19_201b.Fa"))


