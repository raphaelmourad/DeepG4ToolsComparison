# Vincent ROCHER (29/03/2021)
# Convert G4 detector output into tsv format seq score
# required packages
require(Biostrings)
require(tidyverse)

input_table <- snakemake@input[["table"]]
input_fas <- snakemake@input[["fas"]]
output_file <- snakemake@output[[1]]

#read fasta file
res <- mclapply(1:length(input_table),function(i){
  fas = readDNAStringSet(input_fas[[i]])
  
  G4pos.G4detector <- read_delim(input_table[[i]],delim = ",")
  
  tibble(names = names(fas),seq = as.character(fas)) %>% 
    left_join(G4pos.G4detector,by ="seq") %>% dplyr::select(names,score) %>% 
    replace_na(list(score = 0))
},mc.cores=as.integer(snakemake@threads)) %>% bind_rows()

res %>% write_tsv(output_file)