# This code take as input the G4Hunter.R output and get probability to have a G4 based on G4Hunter ranger implementation
# Vincent ROCHER on 06/10/2020

#required packages
library(tidyverse)
library(pROC)

args = commandArgs(trailingOnly=TRUE)

input_table <- args[1]
output_tsv <- args[2]
output_pdf <- str_replace(output_tsv,"tsv","pdf")

#Read dataset
dataset <- read_tsv(input_table)

#Get Labels
labels <- dataset %>% pull(1) %>% str_detect("pos") %>% as.numeric()
score <- dataset %>% pull(2)
#Process ROC curve and AUC
rocRFall <- roc(as.factor( labels),score,ci=T)
aucRF <- pROC::auc(rocRFall)

tibble(file = input_table,AUC = round(aucRF, digits = 3)) %>% write_tsv(output_tsv)

pdf(output_pdf,height = 4, width = 4)
plot(rocRFall, col = "black",main = paste0("AUC = ",
                                            round(aucRF, digits = 3)))
dev.off()