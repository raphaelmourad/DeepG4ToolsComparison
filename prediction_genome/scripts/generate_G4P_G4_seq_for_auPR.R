#Author: Vincent ROCHER
#Date: 24/02/2021
#from GSE133379 QG4 data generate the overlap between G4seq to get the active set (like BG4G4seq)
require(tidyverse)
require(rtracklayer)
require(plyranges)
require(gkmSVM)

PeaksG4seq.GR <- import.bed("/media/ElissarDisk/Elissar/Projects/DeepG4/Data/G4seq/hg19/GSE63874/Data/Bed/Peaks_G4seq_hg19.bed")

G4Pfiles <- "/media/ElissarDisk/Elissar/Projects/DeepG4/Data/GEO/hg19/GSE133379_G4P" %>% 
  list.files(".narrowPeak.gz",full.names = T) 
names(G4Pfiles) <- basename(G4Pfiles) %>% str_remove(".narrowPeak.gz")

G4Pdata <- G4Pfiles %>% map(read_narrowpeaks) %>% 
  map(as_tibble) %>% 
  bind_rows(.id = "exp")

G4Pdata <- G4Pdata %>%
  mutate(exp = str_replace(exp,"DF-1","DF1")) %>% 
  mutate(exp = str_replace(exp,"HeLa-S3","HeLaS3")) %>% 
  separate(exp,into = c("cell_line","experiment","species","replicate"),sep="-")
#Generate BG4G4 only for human data 
G4Pdata.hg19 <- G4Pdata %>% filter(species == "hg19")%>% 
  split(.,.$cell_line)


G4PG4seqdata.hg19 <- G4Pdata.hg19 %>% map(function(x){
  x %>% as_granges()  %>% mutate(start = start+peak)%>% anchor_start() %>% mutate(width = 1) %>% 
    anchor_center() %>% mutate(width = 201) %>% 
    regioneR::filterChromosomes(keep.chr=paste0("chr",c(1:22,"X","Y"))) %>% 
    unique() %>% 
    filter_by_overlaps(PeaksG4seq.GR) 
})

walk(names(G4PG4seqdata.hg19),function(exp){
  write_bed(G4PG4seqdata.hg19[[exp]],str_c("bed/Peaks_G4P_G4seq_",exp,"_hg19_201b.bed"))
})

#Do the same with our own peak calling peaks
PeaksG4seq.GR <- import.bed("/media/ElissarDisk/Elissar/Projects/DeepG4/Data/G4seq/hg19/GSE63874/Data/Bed/Peaks_G4seq_hg19.bed")

G4Pfiles.reforged <- "/media/ElissarDisk/Elissar/Projects/DeepG4/Data/GEO/hg19/GSE133379_G4P/MACS" %>% 
  list.files(".narrowPeak",full.names = T) 
names(G4Pfiles.reforged) <- basename(G4Pfiles.reforged) %>% str_remove(".narrowPeak")

G4Pdata.reforged <- G4Pfiles.reforged %>% map(read_narrowpeaks) %>% 
  map(as_tibble) %>% 
  bind_rows(.id = "exp")

G4Pdata.reforged <- G4Pdata.reforged %>%
  mutate(exp = str_replace(exp,"DF-1","DF1")) %>% 
  mutate(exp = str_replace(exp,"HeLa-S3","HeLaS3")) 
#Generate BG4G4 only for human data 
G4Pdata.hg19 <- G4Pdata.reforged %>% 
  split(.,.$exp)


G4PG4seqdata.hg19 <- G4Pdata.hg19 %>% map(function(x){
  x %>% as_granges()  %>% mutate(start = start+peak)%>% anchor_start() %>% mutate(width = 1) %>% 
    anchor_center() %>% mutate(width = 201) %>% 
    regioneR::filterChromosomes(keep.chr=paste0("chr",c(1:22,"X","Y"))) %>% 
    unique() %>% 
    filter_by_overlaps(PeaksG4seq.GR) 
})

walk(names(G4PG4seqdata.hg19),function(exp){
  write_bed(G4PG4seqdata.hg19[[exp]],str_c("bed/",exp,"_hg19_201b.bed"))
})
