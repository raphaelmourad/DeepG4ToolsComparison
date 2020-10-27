library(Biostrings)
library(DeepG4)
library(rsample)

# Read positive and segative set of sequences 
sequences.pos <- readDNAStringSet(system.file("extdata", "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b.Fa", package = "DeepG4"))
sequences.ctrl <- readDNAStringSet(system.file("extdata", "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM.Fa", package = "DeepG4"))
sequences <- c(sequences.pos,sequences.ctrl)
# Generate classes
Y <- c(rep(1,length(sequences.pos)),rep(0,length(sequences.ctrl)))