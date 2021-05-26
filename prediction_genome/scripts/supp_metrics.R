# Vincent ROCHER (31/03/2021)
# Generate auPR for each dataset
# required packages
if (!require(plyranges))
{
  BiocManager::install("plyranges")
}
require(plyranges)
# if (!require(yardstick))
# {
#   install.packages("yardstick")
# }
require(yardstick)
require(tidyverse)
output_tsv <- snakemake@output[["tsv"]]

readBed <- . %>% read_tsv(col_names = F) %>% dplyr::select(1:3) %>% setNames(c("seqnames","start","end")) %>% as_granges()
bed <- snakemake@input[["bed"]] %>% read_bed()
if(unique(width(bed)) != 201){
  bed <- snakemake@input[["bed"]] %>% readBed
}

dataset <- snakemake@input[["table"]] %>% read_tsv()
bin_genome <- dataset %>% pull(1) %>% GRanges()
ov <- findOverlaps(bin_genome,bed,minoverlap = 50)

labels <- rep(0,length(bin_genome))
labels[queryHits(ov)] <- 1

score <- dataset %>% pull(2)
score <- as.numeric(score > 0.5) %>% as.factor()
datapred <- tibble(labels=factor(labels,levels=c(1,0)),score=factor(score,levels=c(1,0)))
sres <- labels %>% table %>% enframe() %>% spread(key = name,value = value) %>% mutate(percentage_1 = (`1`/`0`)*100)

precision_metrics <- datapred %>% precision(labels,score) %>% pull(.estimate)
accuracy_metrics <- datapred %>% accuracy(labels,score) %>% pull(.estimate)

cbind(sres,tibble(file = snakemake@input[["table"]],FDR = round(1-precision_metrics, digits = 3),
                  accuracy = round(accuracy_metrics,digits = 3))) %>% write_tsv(output_tsv)
