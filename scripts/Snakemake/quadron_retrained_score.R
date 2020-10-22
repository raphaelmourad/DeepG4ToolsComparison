# Generate score for Quadron retrained model on active G4
# Author : Vincent ROCHER
# date : 02/10/20
# From scripts/data_generation/quadron_features_for_quadron_retrained.R
if (!require(xgboost))
{
  install.packages("https://cran.r-project.org/src/contrib/xgboost_1.2.0.1.tar.gz",repos=NULL,type="source")
  if(!require(xgboost)) stop("Package not found")
}
library(xgboost)

library(tidyverse)
library(Biostrings)


args = commandArgs(trailingOnly=TRUE)
model_path <- args[1]
input_file <- args[2]
output_file <- args[3]
Quadron.test <- read_tsv(input_file)

predSet <- xgb.DMatrix(data=as.matrix(Quadron.test %>%  select(-names)))
QModel <- xgb.load(model_path)
Model.Pred <- predict(QModel,predSet , reshape=T)

res <- tibble(names = pull(Quadron.test,names),
              score = Model.Pred)


res %>% write_tsv(output_file)