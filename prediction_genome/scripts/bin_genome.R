# Vincent ROCHER
# created 29/03/2021
# CBI LBCMCP Legube Team, Toulouse
#Generate binned genome data (201)
require(Biostrings)
if (!require(plyranges))
{
  BiocManager::install("plyranges")
}
require(plyranges)
if (!require(BSgenome.Hsapiens.UCSC.hg19))
{
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}
require(BSgenome.Hsapiens.UCSC.hg19)
require(tidyverse)
seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
#GENERATE BIN 201 FROM ALL GENOME

bin_genome <- tileGenome(seqlens[1:23],tilewidth = 201,cut.last.tile.in.chrom = T)
bin_genome <- bin_genome %>% filter(width == 201)


G4pos.seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,bin_genome) 
names(G4pos.seq) <- paste(bin_genome)
resFreq <- Biostrings::letterFrequency(G4pos.seq,"N",as.prob = T)
testNFreq <- as.vector(resFreq>0.1)


bin_genome <- bin_genome[!testNFreq]
G4pos.seq <- G4pos.seq[!testNFreq]

nbSplits <- as.integer(snakemake@params["NBsplit"])
bin_genome_len <- 1:length(bin_genome) %>% as_tibble() %>% mutate(group = ntile(value,nbSplits))



mclapply(1:nbSplits,function(one_split){
  my_indexes <- bin_genome_len %>% filter(group == one_split) %>% pull(value)
  bin_genome[my_indexes] %>% write_bed(str_c(snakemake@params["mydir"],one_split,".bed"))
  writeXStringSet(G4pos.seq[my_indexes],str_c(snakemake@params["mydir"],one_split,".Fa"))
},mc.cores = as.integer(snakemake@threads))
# })



