require(tidyverse)

# files <- list.files("results", recursive=TRUE, full.names=TRUE)
# files <- files[str_detect(files,"AUC/.*\\.tsv")]
res <- snakemake@input  %>% map(read_tsv) %>% map(gather,key = metric,value = value,-file) %>% bind_rows()

fres <- res %>% 
  mutate(Exp = str_extract(file,"GSE[0-9]+|breast-cancer-PDTX")) %>% 
  mutate(size = str_extract(file,"201b|500b|253b")) %>% 
  mutate(Ctrl = str_extract(file,"Ctrl(_[A-Za-z0-9]+|)")) %>% 
  mutate(Tool = str_extract(basename(file),"\\..+")) %>%
  mutate(Tool = str_remove(Tool,".8_42_Ctrl_(gkmSVM|neg).|.")) %>%
  mutate(cell_line = str_extract(file,"qG4|HaCaT|K562|HEKnp|H1975|A549|293T|HeLaS3")) %>% 
  mutate(TypeExp = str_extract(file,"Peaks_G4P_G4seq|TestSet_Peaks_qG4|TestSet_Peaks_G4seqpm_BG4|TestSet_Peaks_BG4_G4seq|Promoters_qG4|Promoters_G4seq_BG4|Promoters_BG4_G4seq|Peaks_G4seqpm_BG4|Peaks_BG4_G4seq|Peaks_G4seq_qG4|Peaks_G4seq_BG4|Peaks_qG4|Peaks_G4seqpm_qG4")) %>% 
  dplyr::select(-file)

#Raf en veut
fres  %>% 
  unite(Ctrl_dat,TypeExp,Exp,cell_line,Ctrl,size) %>% unite(Ctrl_dat,metric,Ctrl_dat) %>% spread(key = Ctrl_dat,value = value) %>% write_tsv(snakemake@output[[1]])
#Raf en veut
fres  %>%  group_by(Exp,Ctrl,cell_line,TypeExp,size,metric) %>% arrange(desc(value)) %>% slice(1) %>% 
  write_tsv(snakemake@output[[2]])

fres %>% unite(Ctrl_dat,TypeExp,Exp,cell_line,Ctrl,size) %>% group_by(Ctrl_dat,metric) %>% arrange(desc(value)) %>% mutate(Classement = 1:dplyr::n()) %>% 
  dplyr::select(-value) %>% spread(key = Ctrl_dat,value = Classement) %>% write_tsv(snakemake@output[[3]])

fres %>% write_tsv(snakemake@output[[4]])