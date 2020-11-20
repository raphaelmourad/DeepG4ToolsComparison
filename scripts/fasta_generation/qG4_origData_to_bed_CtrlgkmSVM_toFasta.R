#' Created on 19 Oct 2020 by Elissar H.N. 
#' CBI, LBCMCP, Legube Team, Toulouse 
#' This code is to generate the bed and fasta for the breastCancer original data and with the overlapping with G4seq 
#' Also for the overlapping between G4seq and the breastCancer data --> so here 3 datasets are generated 

setwd(dir = "/media/ElissarDisk/Elissar/Projects/DeepG4/Pipeline_Data")
pacman::p_load(BSgenome.Hsapiens.UCSC.hg19,TxDb.Hsapiens.UCSC.hg19.knownGene,BSgenome.Hsapiens.UCSC.hg19.masked,
               GenomicRanges, plyranges,readr,tibble,gkmSVM, dplyr,Biostrings,rtracklayer,stringr)


genome <- BSgenome.Hsapiens.UCSC.hg19
Seqnames<- str_c("chr",c(1:22,"X"))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
Proms <- promoters(genes(txdb), upstream = 250, downstream = 250) 
PromsGR <- Proms[(seqnames(Proms) %in% Seqnames)]


Tab <-  ("Orig_Data/breastCancer/qG4-ChIP-seq-of-breast-cancer-PDTX/peak_counts.norm_filtered.tab")
data <- read.table(Tab, header = T, sep = "\t", fill = TRUE)

## Load G4seq Data 
PeaksG4seq.GR <- read.table("Bed/Peaks_G4seq_hg19.bed") %>% setNames(c("seqnames","start","end"))%>% GRanges


 ## Save original peaks length into bed
BedFile <- data %>% dplyr::select(c(chr, start, end)) %>% write_tsv(file = paste0("Bed/Peaks_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19.bed") ,
                                                                    col_names = FALSE )

## reload the bed 
data %>% mutate(diff = end-start) %>% dplyr::select(diff) %>% summary()

### Median is 253 pb , so Start+end/2 +- 126 


Peaks.GR <- data %>% dplyr::select(chr,start,end) %>%  setNames(c("seqnames","start","end")) %>% as_granges() %>% mutate(start = end) %>% anchor_center() %>% mutate(width = 253)

#### Subset with G4seq and write bed 
SBOs<- subsetByOverlaps(Peaks.GR,PeaksG4seq.GR)
SBO <- SBOs %>%  mutate(start = (start+end)/2) %>% mutate(end=start) %>% anchor_center() %>% mutate(width = 201)%>% as.data.frame() %>%  dplyr::select(c("seqnames","start","end")) %>% 
write_tsv(file = paste0("Bed/Peaks_qG4_G4seq_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b.bed"),col_names = FALSE)



## %>% select(c(V1,V2,V3)) %>%  setNames(c("seqnames","start","end")) 
path2 <- paste0("Bed/Peaks_qG4_G4seq_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b.bed") 
PeaksqG4.GR <- path2 %>%  read.table() %>% setNames(c("seqnames","start","end")) %>% as_granges()
Peaks.seq <-getSeq(genome,PeaksqG4.GR)
names(Peaks.seq)=paste0("Peaks_qG4-ChIP-seq_BrC_G4seq_peak",1:length(Peaks.seq))
writeXStringSet(Peaks.seq, paste0("Fasta/Peaks_qG4_G4seq_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b.Fa"))



######## Generate Control gkmSVM 
## open it to use its length later 
input <- read.table(path2)
## generate bed containing matching repeat and GC content as the input bed file(BG4_G4seq)
genNullSeqs(path2,
            genomeVersion = 'BSgenome.Hsapiens.UCSC.hg19.masked',
            outputBedFN= paste0("Bed/Peaks_qG4_G4seq_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM.bed"),
            outputPosFastaFN = 'Tobedeleted.Fa',
            outputNegFastaFN = 'Tobedeleted2.Fa',
            nMaxTrials=20,xfold=16)

## delete the fasta not in need 
file.remove(c("Tobedeleted.Fa","Tobedeleted2.Fa"))

## read the output we generated 
CtrlBed <- read.table( paste0("Bed/Peaks_qG4_G4seq_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM.bed"))
## add 2 into the end position's column to get 201 bp not 199


## transform data into granges 
Unfiltered.Bed.GR <- GRanges(seqnames = CtrlBed$V1,
                             ranges = IRanges(start = CtrlBed$V2,
                                              end = as.numeric(CtrlBed$V3)+1))
## extract only specific chomosomes                                                                                           
Peaks.GR.Filtered <- Unfiltered.Bed.GR[(seqnames(Unfiltered.Bed.GR) %in% Seqnames)]

## write the output by replacing the old one 
write_tsv( as.data.frame(Peaks.GR.Filtered)[1:nrow(input),1:3],
            file = paste0("Bed/Peaks_qG4_G4seq_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM.bed"),
            col_names = FALSE)

## Control gkmSVM to fasta 
CtrlBed <- read.table( paste0("Bed/Peaks_qG4_G4seq_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM.bed")) %>% 
 setNames(c("seqnames", "start","end")) %>% as_granges()
Peaks.seq.Ctrl <-getSeq(genome,CtrlBed)
names(Peaks.seq.Ctrl)=paste0("Peaks_qG4-ChIP-seq_BrC_G4seq_Ctrl_gkmSVM_peak",1:length(Peaks.seq.Ctrl))
writeXStringSet(Peaks.seq.Ctrl, paste0("Fasta/Peaks_qG4_G4seq_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM.Fa"))


##################################################################################################################################
#### OVERLAP G4SEQ_QG4
PeaksG4seq.GR <- read.table("Bed/Peaks_G4seq_hg19_201b.bed") %>% setNames(c("seqname","start","end"))

# PeaksG4seq.GR.253 <- read.table("Data/G4seq/hg19/GSE63874/Data/Bed/Peaks_G4seq_hg19_253b.bed") %>% setNames(c("seqname","start","end"))

# import bed files of BG4 and G4seq  
PeaksqG4.GR <- paste0("Bed/Peaks_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19.bed") %>%
    read_tsv(col_names = F) %>% setNames(c("seqnames","start","end")) %>% GRanges()



# Granges extraction of overlapped sequences between BG4 and G4seq 
SBOs<- subsetByOverlaps((PeaksG4seq.GR %>% GRanges),PeaksqG4.GR) %>% as.data.frame()

# write common peaks into new bed file
write_tsv( SBOs[,1:3],file =paste0("Bed/Peaks_G4seq_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b.bed"),
           col_names = FALSE)
## tO FASTA 

Peaks.seq.Ctrl <-getSeq(genome,(SBOs %>% GRanges))
names(Peaks.seq.Ctrl)=paste0("Peaks_G4seq_qG4-ChIP-seq_BrC_peak",1:length(Peaks.seq.Ctrl))
writeXStringSet(Peaks.seq.Ctrl, paste0("Fasta/Peaks_G4seq_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b.Fa"))



# generate ctrl G4seq-qG4 
inputBedFN <-paste0("Bed/Peaks_G4seq_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b.bed")
## open it to use its length later 
input <- read.table(inputBedFN)
## generate bed containing matching repeat and GC content as the input bed file(BG4_G4seq)
genNullSeqs(inputBedFN,
            genomeVersion = 'BSgenome.Hsapiens.UCSC.hg19.masked',
            outputBedFN= paste0("Bed/Peaks_G4seq_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM.bed"),
            outputPosFastaFN = 'Tobedeleted.Fa',
            outputNegFastaFN = 'Tobedeleted2.Fa',
            nMaxTrials=30,xfold=20)

## delete the fasta not in need 
file.remove(c("Tobedeleted.Fa","Tobedeleted2.Fa"))

## read the output we generated 
CtrlBed <- read.table( paste0("Bed/Peaks_G4seq_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM.bed"))
## add 2 into the end position's column to get 201 bp not 199


## transform data into granges 
Unfiltered.Bed.GR <- GRanges(seqnames = CtrlBed$V1,
                             ranges = IRanges(start = CtrlBed$V2,
                                              end = as.numeric(CtrlBed$V3)+1))
## extract only specific chomosomes                                                                                           
Peaks.GR.Filtered <- Unfiltered.Bed.GR[(seqnames(Unfiltered.Bed.GR) %in% Seqnames)]

## write the output by replacing the old one 
write_tsv(x = as.data.frame(Peaks.GR.Filtered)[1:nrow(input),1:3],
            file = paste0("Bed/Peaks_G4seq_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM.bed"),
            col_names = FALSE)


#### transform to fasta 
CtrlBed <- read.table( paste0("Bed/Peaks_G4seq_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM.bed")) %>% 
  setNames(c("seqnames", "start","end")) %>% as_granges()
Peaks.seq.Ctrl <-getSeq(genome,CtrlBed)
names(Peaks.seq.Ctrl)=paste0("Peaks_G4seq_qG4-ChIP-seq_BrC_201b_Ctrl_gkmSVM_peak",1:length(Peaks.seq.Ctrl))
writeXStringSet(Peaks.seq.Ctrl, paste0("Fasta/Peaks_G4seq_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM.Fa"))
###############   
