# Take as input the csv format concatenation file of qgrs output and generate a tsv with sequences names
# Author : Vincent ROCHER
#date 05/10/20
# Input infos
# Column 1 - ID (order found in sequence).  x.y where x is the primary id, and y is number assigned overlaps.
# For example, all QGRS listed as 2.y overlap QGRS listed with ID 2 - where 2 is the QGRS resulting
# in the highest G-Score in the group.
# Column 2 - Position of the start of the first tetrad (relative to beginning of input sequence)
# Column 3 - Position of the start of the second tetrad (relative to beginning of input sequence)
# Column 4 - Position of the start of the third tetrad (relative to beginning of input sequence)
# Column 5 - Position of the start of the fourth tetrad (relative to beginning of input sequence)
# Column 6 - Number of tetrads
# Column 7 - G-Score
# Column 8 - Sequence
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
#getting option for Quadron algorithm

input_table <- args[1]
ouput_file <- args[2]
type_calc <- args[3]


G4pos.qgrs <- read_csv(input_table,col_names = F)  %>% dplyr::select(X1,X8) %>% group_by(X1) %>% replace_na(list(X8 = 0))
if(type_calc == "mean"){
  G4pos.qgrs <- G4pos.qgrs %>% summarise(score = mean(X8,na.rm = T))
}else if(type_calc == "max"){
  G4pos.qgrs <- G4pos.qgrs %>% summarise(score = max(X8,na.rm = T))
}else{
  G4pos.qgrs <- G4pos.qgrs %>% summarise(score = sum(X8,na.rm = T))
}

G4pos.qgrs %>% dplyr::rename(names = X1) %>% write_tsv(ouput_file)