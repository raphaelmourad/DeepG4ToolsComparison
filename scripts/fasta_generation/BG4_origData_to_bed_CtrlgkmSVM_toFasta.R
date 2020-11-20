# Created by Elissar H.N. and Vincent R.
# Last modification on Tuesday, 21/01/2020 by Elissar H.N.
# CBI LBCMCP Legube Team, Toulouse
# This script helps to generate bed files containing the coordinates of each peak
# with a length of win*2+1 bp starting from summits as input. 


setwd(dir = "/media/ElissarDisk/Elissar/Projects/DeepG4/Pipeline_Data")
pacman::p_load(rtracklayer,plyranges,tidyverse,Biostrings,GenomicRanges,dplyr,tibble,readr,gkmSVM,BSgenome.Hsapiens.UCSC.hg19)


### Choose Cell line and GSE 
GSE <- c("GSE76688","GSE76688","GSE107690","GSE99205")
CL <- c("HaCaT", "HEKnp", "K562","HaCaT")
#### to filtrate bed according to chromosome names 
Seqnames<- str_c("chr",c(1:22,"X"))
genome <- BSgenome.Hsapiens.UCSC.hg19
win = 100
Combination <- c("BG4_G4seq", "G4seq_BG4")




for (i in 1:length(GSE)) {
    # Load BG4 files and merge summits
    lf <- list.files(paste0("Orig_Data/",CL[i],"/",GSE[i]),recursive=T,full.names=T,pattern = "summits.bed")
    
    # map is to apply function to all element of a vector then bid all rows together 
    Summits.BG4= lf %>% map(read_tsv,col_names = F) %>% map(dplyr::select,1:3) %>% bind_rows() %>% 
        setNames(c("seqnames","start","end")) %>% as_granges()
    
    
    #           
    Peaks.BG4 <- Summits.BG4 %>% mutate(start = end) %>% anchor_center() %>% mutate(width = (win*2)+1)
    
    ##  get rid of unwanted chromosomes 
    Peaks.BG4.Filtered <- Peaks.BG4[seqnames(Peaks.BG4) %in% Seqnames]
    
    as_tibble(Peaks.BG4.Filtered)[,1:3] %>% write_tsv(paste0("Bed/Peaks_BG4_",CL[i],"_",
                                                             GSE[i],"_hg19_",win*2+1,"b.bed"),
                                                      col_names = FALSE)
    
    ##### Orig Narrow peak to no 201 base 
    lf <- list.files(paste0("Orig_Data/",CL[i],"/",GSE[i],"/"),recursive=T,full.names=T,pattern = "narrowPeak")
    
    ## read both files, concatenate them, set the column names and get the granges,finally reduce the granges(get rid of comon peaks) 
    Peaks.BG4= lf %>% map(read_tsv,col_names = F) %>% map(dplyr::select,1:3) %>% bind_rows() %>% 
      setNames(c("seqnames","start","end")) %>% as_granges()
    
    ##  get rid of unwanted chromosomes 
    Peaks.BG4.Filtered <- Peaks.BG4[seqnames(Peaks.BG4) %in% Seqnames]
    
    write.table(x = as.data.frame(Peaks.BG4.Filtered)[,1:3],
                file = paste0("Bed/Peaks_BG4_",CL[i],"_",GSE[i],"_hg19.bed") ,
                append = FALSE,
                quote = FALSE, sep = "\t",  eol = "\n", row.names = FALSE, col.names = FALSE)
    
    
    
    
   
    #################### Find overlap BG4_G4seq and write bed 
    G4seq.Peaks <- "Bed/Peaks_G4seq_hg19.bed" %>% read_tsv()%>% setNames(c("seqnames","start","end")) %>% as_granges()
    BG4.peaks <- paste0("Bed/Peaks_BG4_",CL[i],"_",
                        GSE[i],"_hg19_",win*2+1,"b.bed") %>% read_tsv(col_names = FALSE) %>% setNames(c("seqnames","start","end")) %>% as_granges()
    SBOs_BG4G4seq<- subsetByOverlaps(BG4.peaks,G4seq.Peaks) %>% 
    as_tibble(SBOs_BG4G4seq)[,1:3] %>% write_tsv(paste0("Bed/Peaks_BG4_G4seq_",
                                               CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.bed"),col_names = FALSE)
    
    # 1- Form BG4_G4seq 
    Peaks.seq <- getSeq(genome,SBOs_BG4G4seq)
    names(Peaks.seq)=paste0("BG4_G4seq_",CL[i],"_",GSE[i],"_peak",1:length(SBOs_BG4G4seq))
    
    # write fasta files
    writeXStringSet(Peaks.seq, paste0("Fasta/Peaks_BG4_G4seq_",
                                      CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.Fa"))
    
    
    ###########################" G4seq_BG4 overlap #########################"
    PeaksG4seq.GR.201 <- read_tsv("Bed/Peaks_G4seq_hg19_201b.bed") %>% setNames(c("seqnames","start","end"))%>% GRanges
    
    
    # import bed files of BG4 and G4seq  
    PeaksBG4.GR <- paste0("Bed/Peaks_BG4_",CL[i],"_",GSE[i],"_hg19.bed") %>%
        read_tsv(col_names = F) %>% setNames(c("seqnames","start","end")) %>% GRanges
    
    
    # Granges extraction of overlapped sequences between BG4 and G4seq 
    SBOs_G4seqBG4<- subsetByOverlaps(PeaksG4seq.GR.201,PeaksBG4.GR) ## or filter_by_overlaps from plyranges library 
    SBOs_G4seqBG4 <- SBOs_G4seqBG4 %>% as.data.frame() %>% mutate(start = start-1) 
    # write common peaks into new bed file
    write.table(x = as.data.frame(SBOs_G4seqBG4)[,1:3],file =paste0("Bed/Peaks_G4seq_BG4_",
                                                           CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.bed"),
                quote = FALSE, sep = "\t",
                row.names = FALSE,
                col.names = FALSE )
    
    # 1- Form q
    Peaks.seq.G4seqBG4 <- getSeq(genome,(SBOs_G4seqBG4 %>% GRanges()))
    names(Peaks.seq.G4seqBG4)=paste0("G4seq_BG4_",CL[i],"_",GSE[i],"_peak",1:nrow(SBOs_G4seqBG4))
    
    # write fasta files
    writeXStringSet(Peaks.seq.G4seqBG4, paste0("Fasta/Peaks_G4seq_BG4_",
                                      CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.Fa"))
    
    #Generate gkmSVM control

    ## open it to use its length later
    for (Comb  in Combination) {
        inputBedFN <- paste0("Bed/Peaks_",Comb,"_",
                             CL[i],"_",GSE[i],"_hg19_",win*2+1,"b.bed")
        input_BG4_G4seq <- read.table(inputBedFN)
        ## generate bed containing matching repeat and GC content as the input bed file(BG4_G4seq)
        genNullSeqs(inputBedFN,
                    genomeVersion = 'BSgenome.Hsapiens.UCSC.hg19.masked',
                    outputBedFN= paste0("Bed/Peaks_",Comb,"_",
                                        CL[i],"_",GSE[i],"_hg19_",win*2+1,"b_Ctrl_gkmSVM.bed"),
                    outputPosFastaFN = 'Tobedeleted.Fa',
                    outputNegFastaFN = 'Tobedeleted2.Fa',xfold=16,nMaxTrials =30 )

        ## delete the fasta not in need
        file.remove(c("Tobedeleted.Fa","Tobedeleted2.Fa"))

        ## read the output we generated
        CtrlBed <- read.table( paste0("Bed/Peaks_",Comb,"_",
                                      CL[i],"_",GSE[i],"_hg19_",win*2+1,"b_Ctrl_gkmSVM.bed"))


        ## transform data into granges
        Unfiltered.Bed.GR <- GRanges(seqnames = CtrlBed$V1,
                                     ranges = IRanges(start = CtrlBed$V2,
                                                      end = CtrlBed$V3+1))
        ## extract only specific chomosomes
        Peaks.GR.Filtered <- Unfiltered.Bed.GR[(seqnames(Unfiltered.Bed.GR) %in% Seqnames)]

        as_tibble(Peaks.GR.Filtered)[1:nrow(input_BG4_G4seq),1:3]  %>%  write_tsv( paste0("Bed/Peaks_",Comb,"_",
                                                                                          CL[i],"_",GSE[i],"_hg19_",win*2+1,"b_Ctrl_gkmSVM.bed"),col_names = FALSE)


        # Extract fasta from bed coordinates


        # 2- From BG4_G4seq_CtrlgkmSVM
        input_Ctrl <- read.table(paste0("Bed/Peaks_",Comb,"_",
                                        CL[i],"_",GSE[i],"_hg19_",win*2+1,"b_Ctrl_gkmSVM.bed")) %>% setNames(c("seqname","start","end"))
        Peaks.seq <- getSeq(genome,(input_Ctrl %>% GRanges))
        names(Peaks.seq)=paste0("BG4_G4seq_",CL[i],"_",GSE[i],"_peak_Ctrl_gkmSVM",1:nrow(input_Ctrl))

        # write fasta files
        writeXStringSet(Peaks.seq, paste0("Fasta/Peaks_",Comb,"_",
                                          CL[i],"_",GSE[i],"_hg19_",win*2+1,"b_Ctrl_gkmSVM.Fa"))

    }
    
    
    
    
  
    
}

