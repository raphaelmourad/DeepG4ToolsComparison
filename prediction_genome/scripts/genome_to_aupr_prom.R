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

input_genes <- snakemake@params[["genes"]] %>% read_bed() %>% promoters(500,500)

output_tsv <- snakemake@output[["tsv"]]
output_rds <- snakemake@output[["rds"]]
output_pdf <- snakemake@output[["pdf"]]

readBed <- . %>% read_tsv(col_names = F) %>% dplyr::select(1:3) %>% setNames(c("seqnames","start","end")) %>% as_granges()
bed <- snakemake@input[["bed"]] %>% read_bed()
if(unique(width(bed)) != 201){
  bed <- snakemake@input[["bed"]] %>% readBed
}

dataset <- snakemake@input[["table"]] %>% read_tsv()
bin_genome <- dataset %>% pull(1) %>% GRanges()
ov_prom <- findOverlaps(bin_genome,input_genes,minoverlap=100)
bin_genome <- bin_genome[queryHits(ov_prom)]

ov <- findOverlaps(bin_genome,bed,minoverlap = 50)

labels <- rep(0,length(bin_genome))
labels[queryHits(ov)] <- 1

score <- dataset %>% dplyr::slice(queryHits(ov_prom)) %>% pull(2)
datapred <- tibble(labels=factor(labels,levels=c(1,0)),score=score)
sres <- labels %>% table %>% enframe() %>% spread(key = name,value = value) %>% mutate(percentage_1 = (`1`/`0`)*100)

aucRF <- tibble(labels=factor(labels,levels=c(1,0)),score=score) %>% pr_auc(labels,score) %>% pull(.estimate)

data_rds <- datapred %>% pr_curve(labels,score)
data_rds %>% saveRDS(output_rds)
p1 <- data_rds %>% autoplot() +theme_classic() +ggtitle( paste0("AUC = ",round(aucRF, digits = 3)))


cbind(sres,tibble(file = snakemake@input[["table"]],AUC = round(aucRF, digits = 3)) ) %>% write_tsv(output_tsv)
pdf(output_pdf,height = 6, width = 6)
print(p1)
dev.off()