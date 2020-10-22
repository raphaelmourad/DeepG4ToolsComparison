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
##Launch all functions
evaluate_model <- function(y,y.pred.prob,y.pred.classes){
  predict_value <- prediction_table(y,y.pred.prob,y.pred.classes)
  table.1 <- predict_value %>% prediction_metrics
  p.1 <- predict_value %>% plot_roc()
  p.2 <- predict_value %>% plot_confusion_matrix()
  return(list(table.1,p.1,p.2))
}
#Generate coefficients for lasso regression
generateCoeff <- function(cvlasso){
	lambda=cvlasso$lambda.min
	coefLasso=cvlasso$glmnet.fit$beta[,which(cvlasso$lambda==lambda)] %>% enframe() %>%
		mutate(value = abs(value)) %>%
		mutate(name = forcats::fct_reorder(name,value)) 
}
#Plot betas for lasso model
print_betas <- function(cvlasso,top=50){
	coefLasso <- cvlasso %>% generateCoeff() %>% arrange(desc(value))%>% dplyr::slice(1:top)
	coefLasso %>% ggplot(aes(x=name,y=value,fill=value)) + geom_bar(stat="identity") + scale_fill_viridis_c() + coord_flip() + theme_classic() + theme(legend.position = "none") + ggtitle(str_c("Top",top,"Betas for Lasso logistic regression",sep=" ")) + xlab("Betas") + ylab("coefficients")
}
#Plot importance variables for lasso model
print_variable_imp <- function(cvlasso,x_train,top=50){
	coefLasso <- cvlasso %>% generateCoeff() %>%pull(value)
	sds <- apply(x_train, 2, sd)
	std_coefs <- (coefLasso * sds) %>% enframe() %>%
		mutate(value = abs(value)) %>%
		mutate(name = forcats::fct_reorder(name,value)) %>% arrange(desc(value)) %>% dplyr::slice(1:top)
	std_coefs %>% ggplot(aes(x=name,y=value,fill=value)) + geom_bar(stat="identity") + scale_fill_viridis_c() + coord_flip() + theme_classic() + theme(legend.position = "none") + ggtitle(str_c("Top",top,"importance variables for Lasso logistic regression",sep=" ")) + xlab("Variables") + ylab("Importance")
}
#Plot importance variables for ranger
print_variable_imp_rf <- function(RF,top=50){
	std_coefs <- importance(RF) %>% enframe() %>%
		mutate(value = abs(value)) %>%
		mutate(name = forcats::fct_reorder(name,value)) %>% arrange(desc(value)) %>% dplyr::slice(1:top)
	std_coefs %>% ggplot(aes(x=name,y=value,fill=value)) + geom_bar(stat="identity") + scale_fill_viridis_c() + coord_flip() + theme_classic() + theme(legend.position = "none") + ggtitle(str_c("Top",top,"importance variables for Random Forest",sep=" ")) + xlab("Variables") + ylab("Importance")
}

#Given a path, extract and evaluate the model
#Need path for model in hdf5 format, and a .rds file with test data
AUCforModel <- function(model,BG4Data.path="data/Input_for_keras/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_OneHot_train_test.rds"){
  BG4_data <- readRDS(BG4Data.path)
  c(x_train, y_train) %<-% BG4_data$train
  c(x_test, y_test) %<-% BG4_data$test
  vocab_size <- dim(x_train)[3]
  best_model <- load_model_hdf5(model)
  c(pred_prob,pred_class) %<-% get_proba_classes(best_model,x_test,type = "Model")
  
  predict_value <- prediction_table(y_test,pred_prob,pred_class)
  AUC <- predict_value  %>% roc_auc(truth,pred_prob) %>% pull(.estimate) %>% round(3)
  return(list(predict_value,AUC))
}

