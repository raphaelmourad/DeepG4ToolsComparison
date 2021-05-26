# Vincent ROCHER
# created 05/03/2021
# CBI LBCMCP Legube Team, Toulouse
# create X fasta/atac file given a positive bed and control bed

#required packages
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
library(BSgenome.Hsapiens.UCSC.hg19)
require(tidyverse)

binbed <- snakemake@params[["random_regions"]] %>% read_bed()
#Functions
NormBW <- function(x,binbed){
  ranges_bin <- x %>% filter_by_overlaps(binbed) %>% .$score  %>% range(na.rm = TRUE, finite = TRUE)
  x$score <- scales::rescale(x$score,to = c(0,1),from = ranges_bin)
  seqlengths(x) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(x)]
  x <- coverage(x,weight = "score")
  return(x)
}


ATAC_dataset <- snakemake@input[["atac_data"]]
message(ATAC_dataset)
xbw <- ATAC_dataset %>% import.bw %>% NormBW(binbed)
export.bw(xbw,snakemake@output[["atac_normalized"]])