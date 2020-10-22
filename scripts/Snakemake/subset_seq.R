# Take a fasta as input and subset each sequences
# Author Vincent ROCHER
# date : 07/10/20
library(Biostrings)
args = commandArgs(trailingOnly=TRUE)

subseq.size <- as.integer(args[1])
input_fas <- args[2]
output_fas <- args[3]

G4pos.seq <- readDNAStringSet(input_fas)

my_testing <- width(G4pos.seq)>subseq.size

G4pos.seq[my_testing] <- subseq(G4pos.seq[my_testing],1,subseq.size)

writeXStringSet(G4pos.seq,output_fas,width = subseq.size)