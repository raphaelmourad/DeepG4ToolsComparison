# Elissar Nassereddine 
# created on Tuesday, 12/11/2019
# CBI LBCMCP Legube Team, Toulouse
# This script is coded to extract only the overlapping peaks' positions between 
# the g4seq Peaks (plus and minus combined)(in vitro) that overlap with the BG4 peaks(in vivo dataset)

# Set the working Directory 
setwd(dir = "/media/ElissarDisk/Elissar/Projects/DeepG4/Pipeline_Data/")

# Load libraries 
pacman::p_load(BSgenome.Hsapiens.UCSC.hg19,rtracklayer,Biostrings,plyranges,
               tidyverse,dplyr,BSgenome.Hsapiens.UCSC.hg19.masked,gkmSVM,stringr,
               readr,tibble,GenomicRanges,tibble)


Seqnames<- str_c("chr",c(1:22,"X"))
# Load genome in variable 
genome <- BSgenome.Hsapiens.UCSC.hg19

### loop over cell lines 
GSE <- c("GSE76688","GSE76688","GSE107690","GSE99205","qG4-ChIP-seq-of-breast-cancer-PDTX")
CL <- c("HaCaT","HEKnp", "K562","HaCaT","breastCancer")
Comb <- c("BG4","BG4","BG4","BG4","qG4")
win <- 100

## load G4seq plus minus data 
PeaksG4seq.GR <- read.table("Bed/Peaks_G4seqpm_hg19_201b.bed") %>% setNames(c("seqnames","start","end")) %>% GRanges


for (i in 1:length(GSE)){    
  # import bed files of BG4/qG4  
  PeaksBG4.GR <- paste0("Bed/Peaks_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19.bed") %>%
    read_tsv(col_names = F) %>% setNames(c("seqnames","start","end")) %>% GRanges()
  
  
  # Granges extraction of overlapped sequences between G4seq and BG4/qG4
  SBOs<- subsetByOverlaps(PeaksG4seq.GR,PeaksBG4.GR) ## or filter_by_overlaps from plyranges library 

  # Write common peaks into new bed file
  SBOs %>% as_tibble() %>% dplyr::select(seqnames, start, end) %>%  write_tsv(paste0("Bed/Peaks_G4seqpm_",Comb[i],"_",
                                                                                            CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.bed"), col_names = FALSE)

  
  ### write Fastas with headers 
  G4seqpm.seq <- getSeq(genome,(SBOs)) 
  names(G4seqpm.seq) <- paste0(paste0("G4seqpm_BG4_",CL[i],"_",GSE[i],"_Peaks_"),1:nrow(as.data.frame(SBOs)))
  writeXStringSet(G4seqpm.seq, paste0("Fasta/Peaks_G4seqpm_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.Fa"))
  

  
  
  ####  GENERATE THE gkmSVM CONTROL 
  inputBedFN <-paste0("Bed/Peaks_G4seqpm_",Comb[i],"_",
                      CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.bed")
  ## open the posi dataset to use its length later 
  input <- read.table(inputBedFN)
  ## generate bed containing matching repeat and GC content as the input bed file(BG4_G4seq)
  genNullSeqs(inputBedFN,
              genomeVersion = 'BSgenome.Hsapiens.UCSC.hg19.masked',
              outputBedFN= paste0("Bed/Peaks_G4seqpm_",Comb[i],"_",
                                  CL[i],"_",GSE[i],"_hg19_",win*2+1,"b_Ctrl_gkmSVM.bed"),
              outputPosFastaFN = 'Tobedeleted.Fa',
              outputNegFastaFN = 'Tobedeleted2.Fa',xfold=24,nMaxTrials = 34)
  
  ## delete the fasta not in need 
  file.remove(c("Tobedeleted.Fa","Tobedeleted2.Fa"))
  
  ## read the output we generated 
  CtrlBed <- read.table( paste0("Bed/Peaks_G4seqpm_",Comb[i],"_",
                                CL[i],"_",GSE[i],"_hg19_",win*2+1,"b_Ctrl_gkmSVM.bed"))
  
  
  ## transform data into granges 
  Unfiltered.Bed.GR <- GRanges(seqnames = CtrlBed$V1,
                               ranges = IRanges(start = CtrlBed$V2,
                                                end = CtrlBed$V3+1))
  ## extract only specific chomosomes                                                                                           
  Peaks.GR.Filtered <- Unfiltered.Bed.GR[(seqnames(Unfiltered.Bed.GR) %in% Seqnames)]
  
  ## write the output by replacing the old one 
  as_tibble(Peaks.GR.Filtered)[1:nrow(input),1:3] %>%  write_tsv(paste0("Bed/Peaks_G4seqpm_",Comb[i],"_",
                                                                 CL[i],"_",GSE[i],"_hg19_",win*2+1,"b_Ctrl_gkmSVM.bed"),col_names = FALSE)

  
  ## extract control sequences and write in fasta file 
  test <- as.data.frame(Peaks.GR.Filtered)[1:nrow(input),1:3] %>% select(c("seqnames","start","end")) %>% GRanges
  genome <- BSgenome.Hsapiens.UCSC.hg19
  G4seqpm.seq.Ctrl <- getSeq(genome,test) 
  names(G4seqpm.seq.Ctrl) <- paste0(paste0("G4seqpm_BG4_",CL[i],"_",GSE[i],"_Peaks_Ctrl_gkmSVM_"),1:nrow(input))
  writeXStringSet(G4seqpm.seq.Ctrl, paste0("Fasta/Peaks_G4seqpm_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_",win*2+1,"b_Ctrl_gkmSVM.Fa"))
  
}

