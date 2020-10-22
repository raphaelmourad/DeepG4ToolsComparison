# Matthieu Genais
# created on friday, 24/01/2020
# CBI LBCMCP Legube Team, Toulouse
# This code is made launch merging_fasta Snakemake rule
#in order to create X fasta file given a positive fasta and
# X controls

#required packages
require(Biostrings)

#The length's sequence is given into the window parameter
fasta_pos <- readDNAStringSet(snakemake@input[["fasta_pos"]])
fasta_neg <- readDNAStringSet(snakemake@input[["fasta_ctrl"]])

#renaming fasta sequences to keep this information available later
names(fasta_pos) = paste0(names(fasta_pos),"_pos")
names(fasta_neg) = paste0(names(fasta_neg),"_neg")

#merging the to fasta files
fasta_merged=c(fasta_pos,fasta_neg)

#output a Fasta file into Snakemake

writeXStringSet(fasta_merged,snakemake@output[["fasta_merged"]])

