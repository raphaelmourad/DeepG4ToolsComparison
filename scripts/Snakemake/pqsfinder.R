# Matthieu Genais
# created on friday, 24/01/2020
# CBI LBCMCP Legube Team, Toulouse
# This code computes pqsfinder score for different sequences Snakemake rule

#required packages
require(Biostrings)
require(pqsfinder)
require(tidyverse)

#load The merged fasta file
seq <- readDNAStringSet(snakemake@input[["fas"]])

#running pqsFinder on every sequence
pqs <- lapply(names(seq),function(i){
    x <- seq[[i]]
    tibble(names = i,score = as(pqsfinder(x),"GRanges")$score %>% abs %>% max)
}) %>% bind_rows()

#set NA scores to 0 and same for infinite scores
pqs <- pqs %>% mutate(score = ifelse(is.na(score),0,score)) %>% mutate(score = ifelse(is.finite(score),score,0))

#output a tsv file for Snakemake
pqs %>% write_tsv(snakemake@output[["tsv"]])