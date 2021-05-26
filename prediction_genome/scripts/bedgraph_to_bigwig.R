if (!require(plyranges))
{
  BiocManager::install("plyranges")
}
require(plyranges)
if (!require(BSgenome.Hsapiens.UCSC.hg19))
{
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}
library(BSgenome.Hsapiens.UCSC.hg19)
require(tidyverse)
seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
my_bedgraph <- snakemake@input[["bedgraph"]] %>% import.bedGraph()
seqlengths(my_bedgraph) <- seqlens[names(seqlengths(my_bedgraph))]
my_bedgraph %>% write_bigwig(snakemake@output[["bw"]])
