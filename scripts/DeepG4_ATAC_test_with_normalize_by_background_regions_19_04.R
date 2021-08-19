# Author : ROCHER Vincent
# Date: 19/04/2020 (Modified)
#Data from input from ../data_generation/generate_DeepG4_ATAC_test_with_normalize_by_background_regions_19_04.R
library(keras)
library(yardstick)
library(tidyverse)
#Source scripts
#Prediction functions used for evaluate a model
##Return the prediction table based on real values and predicted values
prediction_table <- function(y,y.pred.prob,y.pred.classes){
  tibble(
    truth = as.factor(y) %>% fct_recode(Yes = "1", No = "0"),
    estimate = as.factor(y.pred.classes) %>% fct_recode(Yes = "1", No = "0"),
    pred_prob = y.pred.prob[,1]
  )
}
##Return three metrics : auc, accuracy and kap
prediction_metrics <- function(predict_value){
  AUC <- predict_value  %>% roc_auc(truth,pred_prob)
  some_metrics <- predict_value %>% metrics(truth, estimate)
  rbind(AUC,some_metrics)
}
##Return the roc curve of the model
plot_roc <- function(predict_value){
  predict_value %>% roc_curve(truth,pred_prob) %>% autoplot()
}
##Return the confusion matrix represented as ggplot object
plot_confusion_matrix <- function(predict_value){
  p <- predict_value %>% conf_mat(truth, estimate) %>% .[[1]] %>%
    as_tibble() %>%
    ggplot(aes(Prediction, Truth, fill = n)) +
    geom_tile(show.legend = FALSE) +
    scale_fill_viridis_c() +
    geom_text(aes(label = n), color = "white", alpha = 1, size = 8) +
    labs(
      title = "Confusion matrix"
    ) + theme_minimal(base_size=18)
  p
}
evaluate_model <- function(y,y.pred.prob,y.pred.classes){
  predict_value <- prediction_table(y,y.pred.prob,y.pred.classes)
  table.1 <- predict_value %>% prediction_metrics
  p.1 <- predict_value %>% plot_roc()
  p.2 <- predict_value %>% plot_confusion_matrix()
  return(list(table.1,p.1,p.2))
}
##Obtain the probability
get_proba_classes <- function(model,x_test,type = "Sequential",treshold = 0.5){
  if(type == "Sequential"){
    pred_prob <- model %>% predict_proba(x_test)
    pred_class <- model %>% predict_classes(x_test)
  }else{
    pred_prob <- model %>% predict(x_test)
    pred_class <- ifelse(pred_prob < treshold,0,1)
  }
  return(list(pred_prob,pred_class))
}
# FLAGS
FLAGS <- flags(
  flag_string("activation","relu"),
  flag_numeric("filters1",900),
  flag_numeric("kernel_size1",24),
  flag_numeric("pool_size1",10),
  flag_numeric("dropout1",0),
  flag_numeric("dense_1",100),
  flag_string("optimizer",'rmsprop'),
  flag_numeric("learning_rate",0.001),
  flag_string("loss",'binary_crossentropy'),
  flag_numeric("epoch",20),
  flag_string("input","/home/disleg/Documents/Vincent/DeepG4/scripts/ATACseq_18012020/PIPELINE/DVP_github/DeepG4ToolsComparison/rds/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_rescale_BW_by_bg_5kb_seuil_2")
)
#Load data
OneHotBG4_data <- readRDS(str_c(FLAGS$input,"_OneHot_train_test.rds"))
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
  layer_conv_1d(filters = FLAGS$filters1, kernel_size = FLAGS$kernel_size1, activation = FLAGS$activation,name = "conv_motif") %>%
  layer_average_pooling_1d(pool_size = FLAGS$pool_size1,name ="average_motif_signal") %>%
  layer_global_max_pooling_1d(name = "max_pooling_for_motif")


atac_input <- layer_input(shape = atac_input_shape,name = "atac_input")

main_output <- layer_concatenate(c(motif_output,atac_input)) %>%
  layer_dropout(FLAGS$dropout1) %>%
  layer_dense(FLAGS$dense_1) %>%
  layer_dense(1) %>%
  layer_activation("sigmoid",name = "main_output")
model <- keras_model(
  inputs = c(conv_input,atac_input),
  outputs = c(main_output)
)

if(FLAGS$optimizer == "rmsprop"){
  opt <- optimizer_rmsprop(lr = FLAGS$learning_rate)
}else if(FLAGS$optimizer == "adam"){
  opt <- optimizer_adam(lr = FLAGS$learning_rate)
}else{
  opt <- optimizer_sgd(lr = FLAGS$learning_rate)
}

model %>% compile(
  optimizer = opt,
  loss = FLAGS$loss,
  metrics = c('accuracy')
)

#train the model
history <- model %>% fit(
  x = list(x_train,x_ATAC_train),
  y = y_train,
  epochs = FLAGS$epoch,
  batch_size = 128,
  validation_split = 0.2,
  callbacks = list(
    callback_model_checkpoint("best_model.h5",save_best_only=TRUE)
    # ,callback_early_stopping(restore_best_weights = TRUE,patience = 5)
  )
)
# model %>% save_model_hdf5("model.h5")
#Model evaluation
plot(history)
model <- load_model_hdf5("best_model.h5")
#Model evaluation
plot(history)
pred_prob <- model %>% predict(list(x_test,x_ATAC_test))
pred_class <- ifelse(pred_prob<0.5,0,1)
c(table.1,p.1,p.2) %<-% evaluate_model(y_test,1-pred_prob,pred_class)
table.1
print(p.1)
print(p.2)
