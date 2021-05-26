if (!require(plyranges))
{
  BiocManager::install("plyranges")
}
require(plyranges)
require(tidyverse)

dataset <- snakemake@input[["table"]] %>% read_tsv()
bin_genome <- dataset %>% pull(1) %>% GRanges()
bin_genome$score <- dataset %>% pull(2)

bin_genome %>% write_bed_graph(snakemake@output[["bedgraph"]])