## Created by Elissar Nassereddine on 20, Oct, 2020 
#' CBI , LBCMCP, Legube Team, Toulouse 31400
#' This script generate the Fastas of the promoter regions that overlap with G4 peaks 
#' and the corresponding control which are the proms regions that don't overlap with G4 

pacman::p_load(BSgenome.Hsapiens.UCSC.hg19,TxDb.Hsapiens.UCSC.hg19.knownGene, GenomicRanges, plyranges,readr,tibble, Biostrings)

setwd(dir = "/media/ElissarDisk/Elissar/Projects/DeepG4/Pipeline_Data/")

# To avoid have to type the whole package name every time, we use the variable name txdb
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
Seqnames<- stringr::str_c("chr",c(1:22,"X"))
Proms <- promoters(genes(txdb), upstream = 250, downstream = 250) 
PromsGR <- Proms %>% as.data.frame() %>%  GRanges()
genome <- BSgenome.Hsapiens.UCSC.hg19


GSE <- c("GSE76688","GSE76688","GSE107690","GSE99205","qG4-ChIP-seq-of-breast-cancer-PDTX")
CL <- c("HaCaT","HEKnp", "K562","HaCaT","breastCancer")
Comb <- c("BG4_G4seq","BG4_G4seq","BG4_G4seq","BG4_G4seq","qG4_G4seq")
win <- 100


for (i in 1:length(GSE)){
  bedBG4 <- read.table(paste0("Bed/Peaks_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.bed") )
  
  BG4GR <- bedBG4 %>% setNames(c("chrom","start","end")) %>% GRanges()
  
  RASalt <- subsetByOverlaps(PromsGR,BG4GR,minoverlap = round(win*1.8)) %>%  as.data.frame() %>% dplyr::select(c( "seqnames","start" , "end"))

  RASaltCtrl <- PromsGR[!(PromsGR %over% BG4GR)] %>%  as.data.frame() %>% dplyr::select(c( "seqnames","start" , "end"))
  
  ## Save data bed 
  RASalt %>%  write.table(file =paste0("Bed/Promoters_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_500b.bed"),
                          quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE )
  ## Save Data bed
  RASaltCtrl[1:nrow(RASalt),] %>% write.table(file =paste0("Bed/Promoters_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_500b_Ctrl.bed"),
                             quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE )
  ##
  
  Prom.seq <- getSeq(genome,(RASalt %>% GRanges)) 
  names(Prom.seq) <- paste0("Promoters_",Comb,"_",CL[i],"_",GSE[i],"_peak",1:length(Prom.seq))
  Prom.Ctrl.seq <- getSeq(genome,(RASaltCtrl[1:nrow(RASalt),] %>% GRanges)) 
  names(Prom.Ctrl.seq) <- paste0("Promoters_",Comb,"_",CL[i],"_",GSE[i],"_peak_Ctrl",1:length(Prom.seq))
      
  
  # write fasta files
  writeXStringSet(Prom.seq,paste0("Fasta/Promoters_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_500b.Fa"))
  writeXStringSet(Prom.Ctrl.seq, paste0("Fasta/Promoters_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_500b_Ctrl.Fa"))
}
