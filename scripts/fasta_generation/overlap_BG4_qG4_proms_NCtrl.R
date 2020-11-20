## Created by Elissar Nassereddine on 20, Oct, 2020 
#' CBI , LBCMCP, Legube Team, Toulouse 31400
#' This script generate the Fastas of the G4s that overlap with proms regions 
#' and the corresponding control which are the proms regions that don't overlap with G4 as the same length of G4 peaks (201 or 253 bases )
#' 
  

 
setwd(dir = "/media/ElissarDisk/Elissar/Projects/DeepG4/Pipeline_Data/")
pacman::p_load(BSgenome.Hsapiens.UCSC.hg19,TxDb.Hsapiens.UCSC.hg19.knownGene, GenomicRanges, plyranges,readr,tibble, Biostrings)

# To avoid have to type the whole package name every time, we use the variable name txdb
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
Seqnames<- stringr::str_c("chr",c(1:22,"X"))
Proms <- promoters(genes(txdb), upstream = 250, downstream = 250) 
PromsGR <- Proms[(seqnames(Proms) %in% Seqnames)]
genome <- BSgenome.Hsapiens.UCSC.hg19


## Vectors 
GSE <- c("GSE76688","GSE76688","GSE107690","GSE99205","qG4-ChIP-seq-of-breast-cancer-PDTX")
CL <- c("HaCaT","HEKnp", "K562","HaCaT","breastCancer")
Comb <- c("BG4_G4seq","BG4_G4seq","BG4_G4seq","BG4_G4seq","qG4_G4seq")
win <- 100



for (i in 1:length(GSE)){
    bedBG4 <- read.table(paste0("Bed/Peaks_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.bed") )
    BG4GR <- bedBG4 %>% setNames(c("chrom","start","end")) %>% GRanges()
    
    RASalt <- subsetByOverlaps(BG4GR,PromsGR,minoverlap = win*2+1) %>%  as.data.frame() %>% dplyr::select(c( "seqnames","start" , "end"))
    
    ## Save posi bed 
    RASalt %>% as_tibble() %>% write_tsv(file =paste0("Bed/Promoters_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.bed"),col_names = FALSE)
   
     #### generate control by extracting proms don't contain G4
    PromsNoG4 <- PromsGR[!(PromsGR %over% BG4GR)] %>%  anchor_center()  %>% mutate(width = (win*2)+1)
    Ctrl_proms <- PromsNoG4 %>% as_tibble() %>% dplyr::select(seqnames, start,end) 
    


    ## shufle proms positions 
    rows <- sample(nrow(Ctrl_proms))
    Ctrl_shuff <- Ctrl_proms[rows, ]
    
    ## Extract same data size as posi set 
    Ctrl_inputLength <- Ctrl_shuff[1:nrow(RASalt),]
    
    ## Save ctrl bed 
    Ctrl_inputLength %>% write_tsv(file =paste0("Bed/Promoters_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_",win*2+1,"b_Ctrl.bed") , col_names =FALSE)

  
    
    ### Generate Fastas 
    
    Prom.seq <- getSeq(genome,(RASalt %>% GRanges)) 
    names(Prom.seq) <- paste0("Promoters_",Comb[i],"_",CL[i],"_",GSE[i],"_peak",1:length(Prom.seq))
    
    Prom.Ctrl.seq <- getSeq(genome,(Ctrl_inputLength %>% GRanges)) 
    names(Prom.Ctrl.seq) <- paste0("Promoters_",Comb[i],"_",CL[i],"_",GSE[i],"_peak_Ctrl",1:length(Prom.Ctrl.seq))
    
    
    # write fasta files
    writeXStringSet(Prom.seq,paste0("Fasta/Promoters_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.Fa"))
    
    writeXStringSet(Prom.Ctrl.seq, paste0("Fasta/Promoters_",Comb[i],"_",CL[i],"_",GSE[i],"_hg19_",win*2+1,"b_Ctrl.Fa"))
}

