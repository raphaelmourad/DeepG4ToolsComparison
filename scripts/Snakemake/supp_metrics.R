# This code take as input the G4Hunter.R output and get probability to have a G4 based on G4Hunter ranger implementation
# Vincent ROCHER on 06/10/2020

#required packages
library(tidyverse)
library(yardstick)

args = commandArgs(trailingOnly=TRUE)

input_table <- args[1]
output_tsv <- args[2]

#Read dataset
dataset <- read_tsv(input_table)

#Get Labels
labels <- dataset %>% pull(1) %>% str_detect("pos") %>% as.numeric() %>% as.factor()
score <- dataset %>% pull(2)
score <- ifelse(score>0.1,1,0) %>% as.factor()

precision_metric <- precision_vec(labels,score)
accuracy_metric <- accuracy_vec(labels,score)
tibble(file = input_table,
       FDR = round(1-precision_metric, digits = 3),
       accuracy = round(accuracy_metric,digits = 3)) %>% write_tsv(output_tsv)
