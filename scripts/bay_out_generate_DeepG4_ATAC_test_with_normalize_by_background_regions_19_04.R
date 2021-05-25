require(rBayesianOptimization)
# require(ParBayesianOptimization)
require(keras)
require(tidyverse)
require(tfruns)
setwd("/home/disleg/Documents/Vincent/DeepG4/scripts/ATACseq_18012020/PIPELINE")

RunKeras_DeepG4_ATAC_rescale_BW_sampling <- function(inputfile = "/home/disleg/Documents/Vincent/DeepG4/scripts/ATACseq_18012020/PIPELINE/data/rds/oneHot/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_rescale_BW_by_bg_5kb_seuil_2_OneHot_train_test.rds",filters1=900,kernel_size1=20,pool_size1=10,dropout1=0,dense_1=100,learning_rate=0.001,optimizer = "rmsprop",loss="binary_crossentropy",activation = "relu",epoch=10){
  filename = str_c("callbacks_bay_out_RunKeras_DeepG4_ATAC_rescale_BW_sampling",
                   str_remove(basename(inputfile),".rds"),
                   ".csv")
  message(message(inputfile))
  OneHotBG4_data <- readRDS(inputfile)
  c(x_train,x_ATAC_train, y_train) %<-% OneHotBG4_data$train
  c(x_test,x_ATAC_test, y_test) %<-% OneHotBG4_data$test
  vocab_size <- dim(x_train)[3]
  # x_ATAC_train <- x_ATAC_train[[1]]
  # x_ATAC_test <- x_ATAC_test[[1]]
  # x_ATAC_train <- apply(x_ATAC_train,1,mean)
  # x_ATAC_test <- apply(x_ATAC_test,1,mean)
  
  conv_input_shape <- c(ncol(x_train),vocab_size)
  atac_input_shape <- 1
  
  
  
  #First part of the model : convolution
  conv_input <- layer_input(shape = conv_input_shape,name = "conv_input")
  #Second part of the model : stack of dense taking motif
  motif_output <- conv_input %>%
    layer_conv_1d(filters = filters1, kernel_size = kernel_size1, activation = activation,name = "conv_motif") %>% 
    layer_average_pooling_1d(pool_size = pool_size1,name ="average_motif_signal") %>% 
    layer_global_max_pooling_1d(name = "max_pooling_for_motif")
  
  
  atac_input <- layer_input(shape = atac_input_shape,name = "atac_input")
  
  main_output <- layer_concatenate(c(motif_output,atac_input)) %>% 
    layer_dropout(dropout1) %>% 
    layer_dense(dense_1) %>% 
    layer_dense(1) %>%
    layer_activation("sigmoid",name = "main_output")
  model <- keras_model(
    inputs = c(conv_input,atac_input),
    outputs = c(main_output)
  )
  
  if(optimizer == "rmsprop"){
    opt <- optimizer_rmsprop(lr = learning_rate)
  }else if(optimizer == "adam"){
    opt <- optimizer_adam(lr = learning_rate)
  }else{
    opt <- optimizer_sgd(lr = learning_rate)
  }
  
  model %>% compile(
    optimizer = opt,
    loss = loss,
    metrics = c('accuracy')
  )
  
  #train the model
  history <- model %>% fit(
    x = list(x_train,x_ATAC_train),
    y = y_train,
    epochs = epoch,
    batch_size = 128,
    validation_split = 0.2,
    callbacks = list(
      # callback_model_checkpoint("best_model.h5",save_best_only=TRUE),
      callback_csv_logger(filename)
      # ,callback_early_stopping(restore_best_weights = TRUE,patience = 5)
    )
  )
  res_acc <- read_csv(filename) %>% arrange(desc(val_accuracy)) %>% dplyr::slice(1) %>% pull(val_accuracy)
  return(list(Score = res_acc, Pred = 0))
  
}

bounds_1 <- list(
  filters1 = c(100L,1000L),
  kernel_size1 = c(20L,32L),
  pool_size1 = c(1L,16),
  dropout1 = c(0,0.4),
  dense_1 = c(10L,100L)
)

bay_out_1 <- BayesianOptimization(
  FUN = RunKeras_DeepG4_ATAC_rescale_BW_sampling,
  bounds = bounds_1,
  init_points = 10,
  n_iter = 20,
  acq = "ucb", verbose = T)

bay_out_1$History %>% write_tsv("results/bay_out_RunKeras_DeepG4_ATAC_rescale_BW_by_bg_5kb_seuil_2_19_04_2021_second_run.tsv")

p_baye_res <- read_tsv("results/bay_out_RunKeras_DeepG4_ATAC_rescale_BW_by_bg_5kb_seuil_2_19_04_2021_second_run.tsv")  %>%
  gather(key = "metrics",value = value,-Round,-Value)%>%
  mutate(Best_col = ifelse(Round == 29,"red","black"))%>%
  mutate(Best_size = ifelse(Round == 29,3,1))%>%
  mutate(Best_shape = ifelse(Round == 29,17,16)) %>%
  dplyr::rename(`Accuracy (validation)`="Value") %>%
  mutate(metrics = str_remove(metrics,"1|_1")) %>%
  ggplot(aes(x=value,y=`Accuracy (validation)`,shape=Best_shape,col=Best_col,size=Best_size)) +
  geom_point() +
  geom_smooth(method = "loess") +
  scale_size_identity() + scale_color_identity() +
  facet_wrap(~metrics,nrow=3,scales="free_x") + scale_shape_identity() + xlab("")

require(tfruns)
cc.params <- list(
  filters1 = 500,
  kernel_size1 = 20,
  pool_size1 = 16,
  dropout1 = 0,
  dense_1 = 100
)
runs <- tuning_run("scripts/DeepG4_ATAC_test_with_normalize_by_background_regions_19_04.R",
                   runs_dir = "rds/runs/DeepG4_ATAC_rescale_BW_by_bg_5kb_seuil_2_19_04_2021",
                   flags = append(cc.params,
                                  list("epoch"=20))
)
