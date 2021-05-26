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
my_bw <- snakemake@input[["bw"]] %>% import.bw()
my_bw %>% filter(score >= 0.1) %>%  mutate(score = round(score,2)) %>% write_bed_graph(snakemake@output[["bed"]])
