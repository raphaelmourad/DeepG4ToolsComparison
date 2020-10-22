require(tidyverse)

# files <- list.files("results", recursive=TRUE, full.names=TRUE)
# files <- files[str_detect(files,"AUC/.*\\.tsv")]
res <- snakemake@input  %>% map(read_tsv) %>% bind_rows()

fres <- res %>% 
  mutate(Exp = str_extract(file,"GSE[0-9]+|qG4")) %>% 
  mutate(Ctrl = str_extract(file,"Ctrl(_[A-Za-z0-9]+|)")) %>% 
  mutate(Tool = str_extract(basename(file),"\\..+")) %>%
  mutate(Tool = str_remove(Tool,".8_42_Ctrl_gkmSVM.|.")) %>%
  mutate(cell_line = str_extract(file,"qG4|HaCaT|K562|HEKnp")) %>% 
  mutate(TypeExp = str_extract(file,"TestSet_Peaks_BG4_G4seq|Promoters_G4seq_BG4|Promoters_BG4_G4seq|BG4_G4seq|G4seq_BG4|qG4")) %>% 
  dplyr::select(-file)

#Raf en veut
fres  %>% 
  unite(Ctrl_dat,TypeExp,Exp,cell_line,Ctrl) %>% 
  spread(key = Ctrl_dat,value = AUC) %>% write_tsv(snakemake@output[[1]])
#Raf en veut
fres  %>%  group_by(Exp,Ctrl,cell_line,TypeExp) %>% arrange(desc(AUC)) %>% slice(1) %>% 
  write_tsv(snakemake@output[[2]])

fres %>% unite(Ctrl_dat,TypeExp,Exp,cell_line,Ctrl) %>% group_by(Ctrl_dat) %>% arrange(desc(AUC)) %>% mutate(Classement = 1:dplyr::n()) %>% 
  dplyr::select(-AUC) %>% spread(key = Ctrl_dat,value = Classement) %>% write_tsv(snakemake@output[[3]])

fres %>% write_tsv(snakemake@output[[4]])