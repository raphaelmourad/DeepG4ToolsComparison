# Matthieu Genais
# created on friday, 24/01/2020
# CBI LBCMCP Legube Team, Toulouse
# This code is made to launch Quadron Snakemake rule
# Updated by Vincent ROCHER on 22/09/2020

#required packages
require(tidyverse)
.libPaths(c("lib/Rlib",.libPaths()))
require(foreach)
require(itertools)

# PATH = snakemake@params[["tools"]] %>% as.character()

source(paste0(Quadron_lib,"/lib/PatternFinder.R"))
if(!file.exists(paste0(Quadron_lib,"/Quadron.lib"))){
  message("Building Quadron.lib from source...")
  source("scripts/imports/Quadron/lib/bitcompile.R")
}
load(paste0(Quadron_lib,"/Quadron.lib"))

formatSeq <- function(myseq){
  list(names(myseq),seq = paste(myseq),length = nchar(myseq))
}

#seqlist <- G4pos.seq[1:100]
#test = Quadron2(seqlist = seqlist,score =T)




Quadron2 <- function(seqlist = NULL,
                     nCPU         = 4,
                     parsed.seq   = "",
                     SeqPartitionBy = 1000000,
                     fillwith = 0,
                     score =TRUE){
  ################################################################################
  
  foreach(i=1:length(seqlist)) %dopar% {
    
    seq = formatSeq(seqlist[i])
    
    name=seq[[1]]
    QP <- foreach(j=isplitVector(1:seq$length, chunks=ceiling(seq$length/SeqPartitionBy)),
                  .combine="rbind", .inorder=TRUE) %dopar% {
                    rng <- range(j)
                    # Adding 15nt tolerance range for retrieving complete G4s from split parts:
                    tol.rng <- c( max(1,(rng[1]-15)), min((rng[2]+15),seq$length) )
                    qp <- PatternFinder(seq=substr(seq$seq, start=tol.rng[1], stop=tol.rng[2]),
                                        str.pattern="(([G]{3}[NATGCU]{1,12}){3,}[G]{3})|(([C]{3}[NATGCU]{1,12}){3,}[C]{3})")
                    qp$start.pos <- qp$start.pos+tol.rng[1]-1
                    return(qp)
                  }
    
    
    # will use only:
    # QP$sequence
    # QP$start.pos
    # QP$seq.length
    # QP$strand
    
    ##############################################################################
    QP.num.occ <- length(QP$start.pos)
    
    if(QP.num.occ==0){
      if(score==TRUE){ 
        nRES <- c("names","score")
        RES <- as.data.frame(matrix(c(name,0), ncol=length(nRES), nrow=1))
        names(RES) <- nRES
        return(RES)
      }
      else{
        nRES <- nRES <- c("names",'O.Acont','O.Ccont','O.Gcont','O.Length','OT.triad.AAG','OT.triad.AGA','OT.triad.AGC','OT.triad.AGG','OT.triad.ATA','OT.triad.CAG','OT.triad.CCC','OT.triad.CGG','OT.triad.CTG','OT.triad.GAA','OT.triad.GAG','OT.triad.GCA','OT.triad.GCC','OT.triad.GCG','OT.triad.GCT','OT.triad.GGA','OT.triad.GGC','OT.triad.GGG','OT.triad.GGT','OT.triad.GTG','OT.triad.GTT','OT.triad.TGC','OT.triad.TGG','OT.triad.TGT','OT.triad.TTG','G4.Loop1.efe','G4.Loop2.efe','G4.Loop3.efe','G4.nonG3plus.Loop1.Length','G4.nonG3plus.Loop2.Length','G4.nonG3plus.Loop3.Length','G4.nonG3plus.Nloops','FLANK5.O.Acont','FLANK5.O.Ccont','FLANK5.O.Gcont','FLANK5.OT.triad.AAA','FLANK5.OT.triad.AAG','FLANK5.OT.triad.AAT','FLANK5.OT.triad.ACA','FLANK5.OT.triad.ACC','FLANK5.OT.triad.AGA','FLANK5.OT.triad.AGC','FLANK5.OT.triad.AGG','FLANK5.OT.triad.ATG','FLANK5.OT.triad.CAA','FLANK5.OT.triad.CAC','FLANK5.OT.triad.CAG','FLANK5.OT.triad.CCA','FLANK5.OT.triad.CCC','FLANK5.OT.triad.CCG','FLANK5.OT.triad.CCT','FLANK5.OT.triad.CGC','FLANK5.OT.triad.CGG','FLANK5.OT.triad.CTC','FLANK5.OT.triad.CTG','FLANK5.OT.triad.CTT','FLANK5.OT.triad.GAA','FLANK5.OT.triad.GAG','FLANK5.OT.triad.GCA','FLANK5.OT.triad.GCC','FLANK5.OT.triad.GCG','FLANK5.OT.triad.GCT','FLANK5.OT.triad.GGA','FLANK5.OT.triad.GGC','FLANK5.OT.triad.GGG','FLANK5.OT.triad.GGT','FLANK5.OT.triad.GTG','FLANK5.OT.triad.TCA','FLANK5.OT.triad.TCC','FLANK5.OT.triad.TCT','FLANK5.OT.triad.TGA','FLANK5.OT.triad.TGC','FLANK5.OT.triad.TGG','FLANK5.OT.triad.TGT','FLANK5.OT.triad.TTT','FLANK3.O.Acont','FLANK3.O.Ccont','FLANK3.O.Gcont','FLANK3.OT.triad.AAA','FLANK3.OT.triad.AAG','FLANK3.OT.triad.ACA','FLANK3.OT.triad.AGA','FLANK3.OT.triad.AGC','FLANK3.OT.triad.AGG','FLANK3.OT.triad.AGT','FLANK3.OT.triad.ATA','FLANK3.OT.triad.ATG','FLANK3.OT.triad.CAG','FLANK3.OT.triad.CCA','FLANK3.OT.triad.CCC','FLANK3.OT.triad.CCG','FLANK3.OT.triad.CCT','FLANK3.OT.triad.CGC','FLANK3.OT.triad.CTG','FLANK3.OT.triad.GAA','FLANK3.OT.triad.GAG','FLANK3.OT.triad.GAT','FLANK3.OT.triad.GCA','FLANK3.OT.triad.GCC','FLANK3.OT.triad.GCG','FLANK3.OT.triad.GCT','FLANK3.OT.triad.GGA','FLANK3.OT.triad.GGC','FLANK3.OT.triad.GGG','FLANK3.OT.triad.GGT','FLANK3.OT.triad.GTG','FLANK3.OT.triad.TAC','FLANK3.OT.triad.TAG','FLANK3.OT.triad.TAT','FLANK3.OT.triad.TCT','FLANK3.OT.triad.TGA','FLANK3.OT.triad.TGC','FLANK3.OT.triad.TGG','FLANK3.OT.triad.TGT','FLANK3.OT.triad.TTT')
        RES <- as.data.frame(matrix(fillwith, ncol=length(nRES), nrow=1))
        names(RES) <- nRES
        RES$names = name
        return(RES)
      }
    }
    else{
      # hence there are PQSs found.
      
      #### FOREACH EXECUTION #########
      # Creating the feature holding data frame (RES), which will correspond to each
      # found PQS. In case the flanks are incomplete because of terminal location of
      # G4, the corresponding row will be populated by NAs only.
      
      #RES=foreach(i = 1:pqs.len, .combine="rbind") %dopar% {
      #pqs.len=100  ### *** ### TEST LINE
      RES <- foreach(ch=isplitVector(1:QP.num.occ, chunks=ceiling(QP.num.occ/1000)),
                     .combine="rbind", .inorder=TRUE) %dopar% {
                       
                       count <- 1
                       
                       CHRES <- as.data.frame(matrix(NA, ncol=length(feature.names.for.df), nrow=length(ch)))
                       names(CHRES) <- feature.names.for.df
                       
                       for(i in ch){
                         #*# print(i, quote=F)
                         
                         # Getting the GENOMIC sequences for PQS and flanks. Note, that
                         # here 5' and 3' assignment naming will be correct for only + strand.
                         pqs.seq    <-  QP$sequence[i]
                         flank5.seq <- substr( seq$seq,
                                               start = (QP$start.pos[i]-50),
                                               stop  = (QP$start.pos[i]-1) )
                         flank3.seq <- substr( seq$seq,
                                               start = (QP$start.pos[i]+QP$seq.length[i]),
                                               stop  = (QP$start.pos[i]+QP$seq.length[i]+49) )
                         
                         #-- Continuing only if both flanks are fully present.
                         if(nchar(flank5.seq)==50 & nchar(flank3.seq)==50){
                           
                           # Correcting the sequences for - strand G4s.
                           if(QP$strand[i]=="-"){
                             pqs.seq    <- get.revcomp(pqs.seq)
                             plus.flank5.seq  <- flank5.seq
                             flank5.seq <- get.revcomp(flank3.seq)
                             flank3.seq <- get.revcomp(plus.flank5.seq)
                             rm(plus.flank5.seq)
                           }
                           
                           # Extracting the required features.
                           
                           dfr <- FeatureExtractorQr(seq=pqs.seq, type="PQS", O=TRUE, OT=TRUE, G4=TRUE, NC=FALSE)
                           
                           dfr <- data.frame(dfr)
                           
                           df.f5 <- FeatureExtractorQr(seq=flank5.seq, type="FLANK5", O=TRUE, OT=TRUE, G4=FALSE)
                           df.f5 <- data.frame(df.f5)
                           #*# names(df.f5) <- paste("FLANK5.",names(df.f5),sep="")
                           
                           df.f3 <- FeatureExtractorQr(seq=flank3.seq, type="FLANK3", O=TRUE, OT=TRUE, G4=FALSE)
                           df.f3 <- data.frame(df.f3)
                           #*# names(df.f3) <- paste("FLANK3.",names(df.f3),sep="")
                           
                           dfr <- cbind(dfr, df.f5, df.f3)
                           CHRES[count,] <- dfr
                         }
                         #-- End of "continuing only if both flanks are fully present".
                         
                         count <- count+1
                       }
                       
                       return(CHRES)
                       
                     }
      #### FOREACH EXECUTION DONE ####
      
      
      if(score==TRUE){
        
        # load("./ModelProc.env")
        RES[,feature.names.for.df] <- t( (t(RES[,feature.names.for.df]) - medians[feature.names.for.df]) / sdevs[feature.names.for.df] )
        
        #info <- INFOline(OUT=info, msg=
        #"NOTE: Loading Quadron core...")
        #load("./QuadronML")
        suppressPackageStartupMessages(library(xgboost))
        suppressPackageStartupMessages(library(plyr))
        
        
        PRED <- rep(NA, QP.num.occ)
        non.na.rows <- which(!is.na(RES[,1]))
        if( length(non.na.rows)!=0 ){
          suppressPackageStartupMessages(library(xgboost))
          suppressPackageStartupMessages(library(plyr))
          PRED[non.na.rows] <- predict(QuadronML$finalModel, newdata=xgb.DMatrix(as.matrix(RES[non.na.rows,])))
          nRES <- c("names","score")
          RES <- as.data.frame(matrix(c(name,max(PRED[non.na.rows])), ncol=length(nRES), nrow=1))
          names(RES) <- nRES
          return(RES)
        }
        else{
          nRES <- c("names","score")
          RES <- as.data.frame(matrix(c(name,NA), ncol=length(nRES), nrow=1))
          names(RES) <- nRES
          return(RES)
        }
        
      }
      
      
    } 
    
    
    PRED <- rep(NA, QP.num.occ)
    non.na.rows <- which(!is.na(RES[,1]))
    
    ##############################################################################
    
    if(score==FALSE){
      if(length(non.na.rows)!=0 )
      {
        suppressPackageStartupMessages(library(xgboost))
        suppressPackageStartupMessages(library(plyr))
        PRED[non.na.rows] <- predict(QuadronML$finalModel, newdata=xgb.DMatrix(as.matrix(RES[non.na.rows,])))
        PRED=PRED[non.na.rows]
        RES[,feature.names.for.df] <- t( (t(RES[,feature.names.for.df]) - medians[feature.names.for.df]) / sdevs[feature.names.for.df] )
        RES = RES[non.na.rows,]
        RES = RES[PRED==max(PRED),]
        RES <- cbind("names"=unlist(name),RES)
        return(RES)}
      else{
        nRES <- c("names",'O.Acont','O.Ccont','O.Gcont','O.Length','OT.triad.AAG','OT.triad.AGA','OT.triad.AGC','OT.triad.AGG','OT.triad.ATA','OT.triad.CAG','OT.triad.CCC','OT.triad.CGG','OT.triad.CTG','OT.triad.GAA','OT.triad.GAG','OT.triad.GCA','OT.triad.GCC','OT.triad.GCG','OT.triad.GCT','OT.triad.GGA','OT.triad.GGC','OT.triad.GGG','OT.triad.GGT','OT.triad.GTG','OT.triad.GTT','OT.triad.TGC','OT.triad.TGG','OT.triad.TGT','OT.triad.TTG','G4.Loop1.efe','G4.Loop2.efe','G4.Loop3.efe','G4.nonG3plus.Loop1.Length','G4.nonG3plus.Loop2.Length','G4.nonG3plus.Loop3.Length','G4.nonG3plus.Nloops','FLANK5.O.Acont','FLANK5.O.Ccont','FLANK5.O.Gcont','FLANK5.OT.triad.AAA','FLANK5.OT.triad.AAG','FLANK5.OT.triad.AAT','FLANK5.OT.triad.ACA','FLANK5.OT.triad.ACC','FLANK5.OT.triad.AGA','FLANK5.OT.triad.AGC','FLANK5.OT.triad.AGG','FLANK5.OT.triad.ATG','FLANK5.OT.triad.CAA','FLANK5.OT.triad.CAC','FLANK5.OT.triad.CAG','FLANK5.OT.triad.CCA','FLANK5.OT.triad.CCC','FLANK5.OT.triad.CCG','FLANK5.OT.triad.CCT','FLANK5.OT.triad.CGC','FLANK5.OT.triad.CGG','FLANK5.OT.triad.CTC','FLANK5.OT.triad.CTG','FLANK5.OT.triad.CTT','FLANK5.OT.triad.GAA','FLANK5.OT.triad.GAG','FLANK5.OT.triad.GCA','FLANK5.OT.triad.GCC','FLANK5.OT.triad.GCG','FLANK5.OT.triad.GCT','FLANK5.OT.triad.GGA','FLANK5.OT.triad.GGC','FLANK5.OT.triad.GGG','FLANK5.OT.triad.GGT','FLANK5.OT.triad.GTG','FLANK5.OT.triad.TCA','FLANK5.OT.triad.TCC','FLANK5.OT.triad.TCT','FLANK5.OT.triad.TGA','FLANK5.OT.triad.TGC','FLANK5.OT.triad.TGG','FLANK5.OT.triad.TGT','FLANK5.OT.triad.TTT','FLANK3.O.Acont','FLANK3.O.Ccont','FLANK3.O.Gcont','FLANK3.OT.triad.AAA','FLANK3.OT.triad.AAG','FLANK3.OT.triad.ACA','FLANK3.OT.triad.AGA','FLANK3.OT.triad.AGC','FLANK3.OT.triad.AGG','FLANK3.OT.triad.AGT','FLANK3.OT.triad.ATA','FLANK3.OT.triad.ATG','FLANK3.OT.triad.CAG','FLANK3.OT.triad.CCA','FLANK3.OT.triad.CCC','FLANK3.OT.triad.CCG','FLANK3.OT.triad.CCT','FLANK3.OT.triad.CGC','FLANK3.OT.triad.CTG','FLANK3.OT.triad.GAA','FLANK3.OT.triad.GAG','FLANK3.OT.triad.GAT','FLANK3.OT.triad.GCA','FLANK3.OT.triad.GCC','FLANK3.OT.triad.GCG','FLANK3.OT.triad.GCT','FLANK3.OT.triad.GGA','FLANK3.OT.triad.GGC','FLANK3.OT.triad.GGG','FLANK3.OT.triad.GGT','FLANK3.OT.triad.GTG','FLANK3.OT.triad.TAC','FLANK3.OT.triad.TAG','FLANK3.OT.triad.TAT','FLANK3.OT.triad.TCT','FLANK3.OT.triad.TGA','FLANK3.OT.triad.TGC','FLANK3.OT.triad.TGG','FLANK3.OT.triad.TGT','FLANK3.OT.triad.TTT')
        RES <- as.data.frame(matrix(NA, ncol=length(nRES), nrow=1))
        names(RES) <- nRES
        RES$names=name
        return(RES)
      }}} %>% bind_rows()
  
}

