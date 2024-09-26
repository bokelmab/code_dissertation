##required libraries
library(data.table)
library(dplyr)
library(magrittr)
library(ranger)
library(grf)
library(bartCause)
library(latex2exp)

## required source files
source('code/helper/functions_simulation.R')

## helper functions
calc_ate_het_criteo <- function(p_dt_unif, p_dt_het){
  
  exp_y1_g1 <- mean(p_dt_het[W==1 & var_group == 1,]$Y)
  exp_y0_g1 <- mean(p_dt_het[W==0 & var_group == 1,]$Y)
  exp_y1_g2 <- mean(p_dt_het[W==1 & var_group == 2,]$Y)
  exp_y0_g2 <- mean(p_dt_het[W==0 & var_group == 2,]$Y)
  w1_g2 <- mean(p_dt_unif[W == 1,]$var_group == 2)
  w0_g2 <- mean(p_dt_unif[W == 0,]$var_group == 2)
  exp_y1 <- w1_g2*exp_y1_g2+(1-w1_g2)*exp_y1_g1
  exp_y0 <- w0_g2*exp_y0_g2+(1-w0_g2)*exp_y0_g1
  est_ate_het <- exp_y1-exp_y0
  var_exp_y0_g1<- var(p_dt_het[W==0 & var_group == 1,]$Y)/nrow(p_dt_het[W==0 & var_group == 1,])
  var_exp_y0_g2<- var(p_dt_het[W==0 & var_group == 2,]$Y)/nrow(p_dt_het[W==0 & var_group == 2,])
  var_exp_y1_g1<- var(p_dt_het[W==1 & var_group == 1,]$Y)/nrow(p_dt_het[W==1 & var_group == 1,])
  var_exp_y1_g2<- var(p_dt_het[W==1 & var_group == 2,]$Y)/nrow(p_dt_het[W==1 & var_group == 2,])
  var_ate_het <- (w1_g2^2)*var_exp_y1_g2+((1-w1_g2)^2)*var_exp_y1_g1+(w0_g2^2)*var_exp_y0_g2+((1-w0_g2)^2)*var_exp_y0_g1
  
  return(list(est_ate = est_ate_het, est_var = var_ate_het))
}

calc_ate_het <- function(p_dt_complete, p_dt_het){
  
  pH <- mean(p_dt_complete$var_group == 2)
  
  ate1 <- mean(p_dt_het[W==1 & var_group == 1,]$Y)-mean(p_dt_het[W==0 & var_group == 1,]$Y)
  ate2 <- mean(p_dt_het[W==1 & var_group == 2,]$Y)-mean(p_dt_het[W==0 & var_group == 2,]$Y)
  est_ate_het <- (1-pH)*ate1+pH*ate2
  
  var_ate1 <- var(p_dt_het[W==1 & var_group == 1,]$Y)/nrow(p_dt_het[W==1 & var_group == 1,])+var(p_dt_het[W==0 & var_group == 1,]$Y)/nrow(p_dt_het[W==0 & var_group == 1,])
  var_ate2 <- var(p_dt_het[W==1 & var_group == 2,]$Y)/nrow(p_dt_het[W==1 & var_group == 2,])+var(p_dt_het[W==0 & var_group == 2,]$Y)/nrow(p_dt_het[W==0 & var_group == 2,])
  var_ate_het <- ((1-pH)^2)*var_ate1+(pH^2)*var_ate2
  
  return(list(est_ate = est_ate_het, est_var = var_ate_het))
}


## set seed
set.seed(121122)

## choose target
target <- 'conversion'
prob_W <- 0.85

## read Criteo data
data <- fread('data/criteo-uplift-v2.1.csv')

## delete menth closing campaign and encode treatment
names(data)[which(names(data)== 'treatment')] <- 'W'
y_data <- data[, target, with = F]
data$Y <- data[, target, with = F]
data <- data[,-c('conversion', 'visit', 'exposure'), with = F]
data$Y <- y_data

## split in training and test set
idx_base <- sample(1:nrow(data), round(0.01*nrow(data)), replace = F)
dt_base <- data[idx_base,]
dt_exp <- data[-idx_base,]

## build random forest models 
rf_base <- ranger(Y ~ ., data = dt_base[W==0,-c('W'), with = F], max.depth = 5, num.trees = 500)
rf1_base <- ranger(Y ~ ., data = dt_base[W==1,-c('W'), with = F], max.depth = 5, num.trees = 500)
rf0_base <- ranger(Y ~ ., data = dt_base[W==0,-c('W'), with = F], max.depth = 5, num.trees = 500)

## obtain sampling parameters
check <- copy(dt_exp[sample(1:nrow(dt_exp), 100000),])
check$pred_prob <- predict(rf_base, check)$predictions

opt_split_data1 <- find_split(check)
plot((0:100)/100,1-opt_split_data1$var_rel_opt, type = 'l', ylim = c(0,1-min(opt_split_data1$var_rel_opt)), ylab = '% variance reduction', 
     xlab = TeX('$p_L$'), lwd = 3, main = 'Criteo')
lines((0:100)/100, 1-opt_split_data1$var_rel_prop, col = 'red', lwd = 3)

pL <- 5/4*(which.min(opt_split_data1$var_rel_opt)-1)/100-1/4

pred_var1 <- mean(check[pred_prob <= quantile(pred_prob, pL),]$pred_prob)*(1-mean(check[pred_prob <= quantile(pred_prob, pL),]$pred_prob))
pred_var2 <- mean(check[pred_prob > quantile(pred_prob, pL),]$pred_prob)*(1-mean(check[pred_prob > quantile(pred_prob, pL),]$pred_prob))

share_high <- (1-pL)*sqrt(pred_var2)/((1-pL)*sqrt(pred_var2)+pL*sqrt(pred_var1))
share_high <- 3/4*share_high+1/4*(1-pL)
#saveRDS(list(pL = pL, share_high = share_high), paste0('results/', 'criteo', '/sampling_parameters.RDS'))

## define strata in experimental data
pred_prob <- c()
for(i_fold in 1:13){
  pred_prob <- c(pred_prob, predict(rf_base, dt_exp[(1:1000000)+(i_fold-1)*1000000,])$predictions)
  print(i_fold)
}
pred_prob <- c(pred_prob, predict(rf_base, dt_exp[13000001:nrow(dt_exp),])$predictions)
dt_exp$var_group <- ifelse(pred_prob > quantile(pred_prob,pL), 2, 1)

## training of uplift model
idx_train <- sample(1:nrow(dt_exp), 100000, replace = F)
dt_train <- dt_exp[idx_train,]
dt_test <- dt_exp[-idx_train,]

rf0 <- ranger(Y ~ ., data = dt_train[W==0,-c('W','var_group'), with = F], max.depth = 5)
rf1 <- ranger(Y ~ ., data = dt_train[W==1,-c('W','var_group'), with = F], max.depth = 5)

## make predictions on test set
pred_t <- c()
phi_hat <- c()
for(i_fold in 1:13){
  pred_t <- c(pred_t, predict(rf1, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions - predict(rf0, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions)
  phi_hat <- c(phi_hat, predict(rf0_base, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions+(1-prob_W)*predict(rf1_base, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions) ## add phihat predictions for covariate adjustment to experimental data
  print(i_fold)
}
pred_t <- c(pred_t, predict(rf1, dt_test[13000001:nrow(dt_test),])$predictions - predict(rf0, dt_test[13000001:nrow(dt_test),])$predictions)
phi_hat <- c(phi_hat, predict(rf0_base, dt_test[13000001:nrow(dt_test),])$predictions+(1-prob_W)*predict(rf1_base, dt_test[13000001:nrow(dt_test),])$predictions) ## add phihat predictions for covariate adjustment to experimental data
dt_test$predt <- pred_t
dt_test$phi_hat <- phi_hat

## set Y to double
dt_test$Y %<>% as.numeric

## repeat sampling process 1000 times
qini_table_het <- NULL
qini_table_unif <- NULL
qini_table_ca <- NULL
qini_table_hsca <- NULL
var_qini_table_het <- NULL
var_qini_table_unif <- NULL
var_qini_table_ca <- NULL
var_qini_table_hsca <- NULL
est_ate_unif <- c()
est_ate_het <- c()
est_ate_ca <- c()
est_ate_hsca <- c()
var_ate_unif <- c()
var_ate_het <- c()
var_ate_ca <- c()
var_ate_hsca <- c()

for(i_boot in 1:1000){
  
  ## obtain HS sample and random sample
  idx_var_group2 <- which(dt_test$var_group == 2)
  idx_var_group1 <- sample(which(dt_test$var_group == 1), (1/share_high-1)*length(idx_var_group2))
  
  test_unif <- dt_test[sample(1:nrow(dt_test), length(idx_var_group1)+length(idx_var_group2)),]
  test_het <- dt_test[c(idx_var_group1, idx_var_group2),]
  
  ## Obtain ATE estimates
  est_ate_unif <- c(est_ate_unif, mean(test_unif[W==1,]$Y)-mean(test_unif[W==0,]$Y))
  est_ate_ca <- c(est_ate_ca, mean(test_unif[W==1,]$Y-test_unif[W==1,]$phi_hat)-mean(test_unif[W==0,]$Y-test_unif[W==0,]$phi_hat))
  var_ate_unif <- c(var_ate_unif, var(test_unif[W==1,]$Y)/nrow(test_unif[W==1,])+var(test_unif[W==0,]$Y)/nrow(test_unif[W==0,]))
  var_ate_ca <- c(var_ate_ca, var(test_unif[W==1,]$Y-test_unif[W==1,]$phi_hat)/nrow(test_unif[W==1,])+var(test_unif[W==0,]$Y-test_unif[W==0,]$phi_hat)/nrow(test_unif[W==0,]))
  
  ## ATE estimate for HS sampling
  res_het <- calc_ate_het(p_dt_complete = dt_test, p_dt_het = test_het)
  est_ate_het <- c(est_ate_het, res_het$est_ate)
  var_ate_het <- c(var_ate_het, res_het$est_var)
  
  ## ATE estimate for HSCA
  res_hsca <- calc_ate_het(p_dt_complete = dt_test, p_dt_het = copy(test_het)[, Y:= Y-phi_hat])
  est_ate_hsca <- c(est_ate_hsca, res_hsca$est_ate)
  var_ate_hsca <- c(var_ate_hsca, res_hsca$est_var)
  
  ## ATE treatment prio
  qini_unif <- c()
  qini_het <- c()
  qini_ca <- c()
  qini_hsca <- c()
  var_qini_unif <- c()
  var_qini_het <- c()
  var_qini_ca <- c()
  var_qini_hsca <- c()
  for(i_cut in 9:0){
    
    cutoff <- i_cut/10
    quant <- quantile(dt_test$predt, cutoff)
    
    ## ate estimation for t-learner on uniform data
    ate_unif <- mean(test_unif[predt >= quant & W == 1,]$Y)-mean(test_unif[predt >= quant & W == 0,]$Y)
    var_unif <- var(test_unif[predt >= quant & W == 1,]$Y)/nrow(test_unif[predt >= quant & W == 1,])+var(test_unif[predt >= quant & W == 0,]$Y)/nrow(test_unif[predt >= quant & W == 0,])
    ate_ca <- mean(test_unif[predt >= quant & W==1,]$Y-test_unif[predt >= quant & W==1,]$phi_hat)-mean(test_unif[predt >= quant & W==0,]$Y-test_unif[predt >= quant & W==0,]$phi_hat)
    var_ca <- var(test_unif[predt >= quant & W==1,]$Y-test_unif[predt >= quant & W==1,]$phi_hat)/nrow(test_unif[predt >= quant & W==1,])+var(test_unif[predt >= quant & W==0,]$Y-test_unif[predt >= quant & W==0,]$phi_hat)/nrow(test_unif[predt >= quant & W==0,])
    res_het <- calc_ate_het(p_dt_complete = dt_test[predt >= quant,], p_dt_het = test_het[predt >= quant,])
    ate_het <- res_het$est_ate
    var_het <- res_het$est_var
    test_het_quant <- copy(test_het[predt >= quant,])
    test_het_quant[, Y := Y-phi_hat]
    res_hsca <- calc_ate_het(p_dt_complete = dt_test[predt >= quant,], p_dt_het = test_het_quant)
    ate_hsca <- res_hsca$est_ate
    var_hsca <- res_hsca$est_var
    
    qini_unif <- c(qini_unif, ate_unif*(10-i_cut))
    qini_het <- c(qini_het, ate_het*(10-i_cut))
    qini_ca <- c(qini_ca, ate_ca*(10-i_cut))
    qini_hsca <- c(qini_hsca, ate_hsca*(10-i_cut))
    var_qini_unif <- c(var_qini_unif, var_unif/100*(10-i_cut)^2)
    var_qini_het <- c(var_qini_het, var_het/100*(10-i_cut)^2)
    var_qini_ca <- c(var_qini_ca, var_ca/100*(10-i_cut)^2)
    var_qini_hsca <- c(var_qini_hsca, var_hsca/100*(10-i_cut)^2)
    
  }
  qini_table_unif %<>% rbind(qini_unif)
  qini_table_het %<>% rbind(qini_het)
  qini_table_ca %<>% rbind(qini_ca)
  qini_table_hsca %<>% rbind(qini_hsca)
  var_qini_table_unif %<>% rbind(var_qini_unif)
  var_qini_table_het %<>% rbind(var_qini_het)
  var_qini_table_ca %<>% rbind(var_qini_ca)
  var_qini_table_hsca %<>% rbind(var_qini_hsca)
  print(i_boot)
  
  ## save results
  # saveRDS(est_ate_ca, paste0('results/', 'criteo', '/est_ate_ca.RDS'))
  # saveRDS(est_ate_het, paste0('results/', 'criteo', '/est_ate_het.RDS'))
  # saveRDS(est_ate_hsca, paste0('results/', 'criteo', '/est_ate_hsca.RDS'))
  # saveRDS(est_ate_unif, paste0('results/', 'criteo', '/est_ate_unif.RDS'))
  # 
  # saveRDS(var_ate_ca, paste0('results/', 'criteo', '/var_ate_ca.RDS'))
  # saveRDS(var_ate_het, paste0('results/', 'criteo', '/var_ate_het.RDS'))
  # saveRDS(var_ate_hsca, paste0('results/', 'criteo', '/var_ate_hsca.RDS'))
  # saveRDS(var_ate_unif, paste0('results/', 'criteo', '/var_ate_unif.RDS'))
  # 
  # saveRDS(qini_table_ca, paste0('results/', 'criteo', '/qini_table_ca.RDS'))
  # saveRDS(qini_table_het, paste0('results/', 'criteo', '/qini_table_het.RDS'))
  # saveRDS(qini_table_hsca, paste0('results/', 'criteo', '/qini_table_hsca.RDS'))
  # saveRDS(qini_table_unif, paste0('results/', 'criteo', '/qini_table_unif.RDS'))
  # 
  # saveRDS(var_qini_table_ca, paste0('results/', 'criteo', '/var_qini_table_ca.RDS'))
  # saveRDS(var_qini_table_het, paste0('results/', 'criteo', '/var_qini_table_het.RDS'))
  # saveRDS(var_qini_table_hsca, paste0('results/', 'criteo', '/var_qini_table_hsca.RDS'))
  # saveRDS(var_qini_table_unif, paste0('results/', 'criteo', '/var_qini_table_unif.RDS'))
  
  print(round((1-mean(var_ate_ca/var_ate_unif))*100,2))
  print(round((1-mean(var_ate_het/var_ate_unif))*100,2))
  print(round((1-mean(var_ate_hsca/var_ate_unif))*100,2))
}

round((1-apply(var_qini_table_ca,2,mean)/apply(var_qini_table_unif,2,mean))*100,2)
round((1-apply(var_qini_table_het,2,mean)/apply(var_qini_table_unif,2,mean))*100,2)
round((1-apply(var_qini_table_hsca,2,mean)/apply(var_qini_table_unif,2,mean))*100,2)

round((1-mean(var_ate_ca/var_ate_unif))*100,2)
round((1-mean(var_ate_het/var_ate_unif))*100,2)
round((1-mean(var_ate_hsca/var_ate_unif))*100,2)


############### Training part ##################################################

calc_auq <- function(p_preds, p_dt_test){
  ecdf_preds <- ecdf(p_preds)
  p_dt_test$dec <- ceiling(ecdf_preds(p_preds)*10)
  qini <- c()
  for(i_dec in 10:1){
    qini <- c(qini, (mean(p_dt_test[dec >= i_dec & W== 1,]$Y)-mean(p_dt_test[dec >= i_dec & W== 0,]$Y))*(11-i_dec)/10)
  }
  
  return(mean(qini))
}


auqct_unif_rf <- c()
auqcs_unif_rf <- c()
auqcx_unif_rf <- c()
auqct_het_rf <- c()
auqcs_het_rf <- c()
auqcx_het_rf <- c()
auqccf_het <- c()
auqccf_unif <- c()
auqccb_het <- c()
auqccb_unif <- c()
for(i_sim in 1:1000){
  
  idx_unif <- sample(1:nrow(dt_exp), 20000)
  idx_var_group1 <- sample(intersect(which(dt_exp$var_group == 1), idx_unif), round((1-share_high)*20000))
  idx_var_group2 <- union(intersect(which(dt_exp$var_group == 2), idx_unif) , sample(setdiff(which(dt_exp$var_group == 2),idx_unif), 20000*share_high-length(intersect(which(dt_exp$var_group == 2), idx_unif))))
  idx_het <- union(idx_var_group1, idx_var_group2)
  train_unif <- dt_exp[idx_unif,]
  train_het <- dt_exp[idx_het,]
  idx_test <- setdiff(which(dt_exp$var_group==2), idx_var_group2)
  idx_test %<>% union(sample(setdiff(which(dt_exp$var_group==1), idx_unif), pL/(1-pL)*length(idx_test)))
  dt_test <- dt_exp[idx_test,]
  dt_test <- dt_test[sample(1:nrow(dt_test), 1000000),] ## cut run time
  
  
  ############ Training on uniformly sampled data ##############################

  ## T-learner
  rf0_unif <- ranger(Y ~ ., data = train_unif[W==0,-c('W', 'var_group'), with = F], max.depth = 5)
  rf1_unif <- ranger(Y ~ ., data = train_unif[W==1,-c('W', 'var_group'), with = F], max.depth = 5)
  predt_unif <- predict(rf1_unif, data = dt_test)$predictions-predict(rf0_unif, data = dt_test)$predictions
  auqct_unif_rf <- c(auqct_unif_rf, calc_auq(predt_unif, dt_test))
  
  ## S-Learner
  rfs_unif <- ranger(Y ~ ., data = train_unif[,-c('var_group'), with = F], max.depth = 5)
  preds_unif <- predict(rfs_unif, data = copy(dt_test)[,W:=1])$predictions-predict(rfs_unif, data = copy(dt_test)[,W:=0])$predictions
  auqcs_unif_rf <- c(auqcs_unif_rf, calc_auq(preds_unif, dt_test))
  
  ## X-Learner
  train_unif_x <- copy(train_unif[,-c('var_group'), with = F])
  train_unif_x1 <- train_unif_x[W==1,]
  train_unif_x0 <- train_unif_x[W==0,]
  train_unif_x1$Y <- train_unif_x1$Y - predict(rf0_unif, train_unif_x1[,-c('W'),with = F])$predictions
  train_unif_x0$Y <- predict(rf1_unif, train_unif_x0[,-c('W'),with = F])$predictions - train_unif_x0$Y
  rfX0_unif <- ranger(Y ~ ., data = train_unif_x0[,-c('W'), with = F], max.depth = 5)
  rfX1_unif <- ranger(Y ~ ., data = train_unif_x1[,-c('W'), with = F], max.depth = 5)
  predx_unif <- 0.15*predict(rfX0_unif, data = dt_test)$predictions+0.85*predict(rfX1_unif, data = dt_test)$predictions
  auqcx_unif_rf <- c(auqcx_unif_rf, calc_auq(predx_unif, dt_test))

  ## Causal random forest
  cf_unif <- causal_forest(X = train_unif[,-c('W', 'Y','var_group'), with = F], Y = train_unif$Y, W = train_unif$W, W.hat = 0.5, num.trees = 500)
  predcf_unif <- predict(cf_unif, dt_test[,-c('W', 'Y','var_group'), with = F])$predictions
  ecdfcf_unif <- ecdf(predcf_unif)
  auqccf_unif <- c(auqccf_unif, calc_auq(predcf_unif, dt_test))
  
  ## causal bart
  cbart <- bartc(response = train_unif$Y, treatment = train_unif$W, 
                 confounders = train_unif[,-c('W', 'Y','var_group'), with = F], 
                 n.samples = 100L, n.burn = 15L, n.chains = 2L, method.trt = rep(0.5, nrow(train_unif)),keepTrees = T)
  predccb_unif <- apply(predict(cbart, newdata = copy(dt_test)[, ps:= 1], type = 'icate'),2,mean)
  ecdfcb_unif <- ecdf(predccb_unif)
  auqccb_unif <- c(auqccb_unif, calc_auq(predccb_unif, dt_test))
  
  ########## Training on heteroskedasticity aware sampled data ##################
  
  ## T-learner
  rf0_het <- ranger(Y ~ ., data = train_het[W==0,-c('W', 'var_group'), with = F], max.depth = 5)
  rf1_het <- ranger(Y ~ ., data = train_het[W==1,-c('W', 'var_group'), with = F], max.depth = 5)
  predt_het <- predict(rf1_het, data = dt_test)$predictions-predict(rf0_het, data = dt_test)$predictions
  auqct_het_rf <- c(auqct_het_rf, calc_auq(predt_het, dt_test))
  
  ## S-Learner
  rfs_het <- ranger(Y ~ ., data = train_het[,-c('var_group'), with = F], max.depth = 5)
  preds_het <- predict(rfs_het, data = copy(dt_test)[,W:=1])$predictions-predict(rfs_het, data = copy(dt_test)[,W:=0])$predictions
  auqcs_het_rf <- c(auqcs_het_rf, calc_auq(preds_het, dt_test))
  
  ## X-Learner
  train_het_x <- copy(train_het[,-c('var_group'), with = F])
  train_het_x1 <- train_het_x[W==1,]
  train_het_x0 <- train_het_x[W==0,]
  train_het_x1$Y <- train_het_x1$Y - predict(rf0_het, train_het_x1[,-c('W'),with = F])$predictions
  train_het_x0$Y <- predict(rf1_het, train_het_x0[,-c('W'),with = F])$predictions - train_het_x0$Y
  rfX0_het <- ranger(Y ~ ., data = train_het_x0[,-c('W'), with = F], max.depth = 5)
  rfX1_het <- ranger(Y ~ ., data = train_het_x1[,-c('W'), with = F], max.depth = 5)
  predx_het <- 0.15*predict(rfX0_het, data = dt_test)$predictions+0.85*predict(rfX1_het, data = dt_test)$predictions
  auqcx_het_rf <- c(auqcx_het_rf, calc_auq(predx_het, dt_test))

  ## Causal random forest
  cf_het <- causal_forest(X = train_het[,-c('W', 'Y','var_group'), with = F], Y = train_het$Y, W = train_het$W, W.hat = 0.5, num.trees = 500)
  predcf_het <- predict(cf_het, dt_test[,-c('W', 'Y','var_group'), with = F])$predictions
  ecdfcf_het <- ecdf(predcf_het)
  auqccf_het <- c(auqccf_het, calc_auq(predcf_het, dt_test))
  
  ## Causal Bart
  cbart_het <- bartc(response = train_het$Y, treatment = train_het$W, 
                     confounders = train_het[,-c('W', 'Y','var_group'), with = F], 
                     n.samples = 100L, n.burn = 15L, n.chains = 2L, method.trt = rep(0.5, nrow(train_het)),keepTrees = T)
  predccb_het <- apply(predict(cbart_het, newdata = copy(dt_test)[, ps:= 1], type = 'icate'),2,mean)
  ecdfcb_het <- ecdf(predccb_het)
  auqccb_het <- c(auqccb_het, calc_auq(predccb_het, dt_test))
  
  print(i_sim)
  
  ## evaluation T-learner rf
  100*(mean(auqct_het_rf/auqct_unif_rf)-1)
  print(paste0('[', round(mean(auqct_het_rf/auqct_unif_rf)-1-1.96*sd(auqct_het_rf/auqct_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100, ';', round(mean(auqct_het_rf/auqct_unif_rf)-1+1.96*sd(auqct_het_rf/auqct_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100,']'))

  ## evaluation S-learner rf
  100*(mean(auqcs_het_rf/auqcs_unif_rf)-1)
  print(paste0('[', round(mean(auqcs_het_rf/auqcs_unif_rf)-1-1.96*sd(auqcs_het_rf/auqcs_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100, ';', round(mean(auqcs_het_rf/auqcs_unif_rf)-1+1.96*sd(auqcs_het_rf/auqcs_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100,']'))

  ## evaluation X-learner rf
  100*(mean(auqcx_het_rf/auqcx_unif_rf)-1)
  print(paste0('[', round(mean(auqcx_het_rf/auqcx_unif_rf)-1-1.96*sd(auqcx_het_rf/auqcx_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100, ';', round(mean(auqcx_het_rf/auqcx_unif_rf)-1+1.96*sd(auqcx_het_rf/auqcx_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100,']'))

  ## evaluation causal random forest
  100*(mean(auqccf_het/auqccf_unif)-1)
  print(paste0('[', round(mean(auqccf_het/auqccf_unif)-1-1.96*sd(auqccf_het/auqccf_unif)/sqrt(length(auqccf_het)),4)*100, ';', round(mean(auqccf_het/auqccf_unif)-1+1.96*sd(auqccf_het/auqccf_unif)/sqrt(length(auqccf_het)),4)*100,']'))

  ## evaluation causal bart
  100*(mean(auqccb_het/auqccb_unif)-1)
  print(paste0('[', round(mean(auqccb_het/auqccb_unif)-1-1.96*sd(auqccb_het/auqccb_unif)/sqrt(length(auqccb_het)),4)*100, ';', round(mean(auqccb_het/auqccb_unif)-1+1.96*sd(auqccb_het/auqccb_unif)/sqrt(length(auqccb_het)),4)*100,']'))
  
  print(mean(auqct_het_rf))
  print(mean(auqcs_het_rf))
  print(mean(auqcx_het_rf))
  print(mean(auqccf_het))
  print(mean(auqccb_het))
  
  ## save results
  # saveRDS(auqct_het_rf, paste0('results/', 'criteo', '/auqct_het_rf.RDS'))
  # saveRDS(auqcs_het_rf, paste0('results/', 'criteo', '/auqcs_het_rf.RDS'))
  # saveRDS(auqcx_het_rf, paste0('results/', 'criteo', '/auqcx_het_rf.RDS'))
  # saveRDS(auqccf_het, paste0('results/', 'criteo', '/auqccf_het.RDS'))
  # saveRDS(auqccb_het, paste0('results/', 'criteo', '/auqccb_het.RDS'))
  
  # saveRDS(auqct_unif_rf, paste0('results/', 'criteo', '/auqct_unif_rf.RDS'))
  # saveRDS(auqcs_unif_rf, paste0('results/', 'criteo', '/auqcs_unif_rf.RDS'))
  # saveRDS(auqcx_unif_rf, paste0('results/', 'criteo', '/auqcx_unif_rf.RDS'))
  # saveRDS(auqccf_unif, paste0('results/', 'criteo', '/auqccf_unif.RDS'))
  # saveRDS(auqccb_unif, paste0('results/', 'criteo', '/auqccb_unif.RDS'))
  
}

#fwrite(dt_exp, paste0('results/', 'criteo', '/dt_exp.csv'))
#fwrite(dt_base, paste0('results/', 'criteo', '/dt_base.csv'))

## calculate results for paper
round(100*(mean(auqct_het_rf/auqct_unif_rf)-1),2)
print(paste0('[', round(mean(auqct_het_rf/auqct_unif_rf)-1-1.96*sd(auqct_het_rf/auqct_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100, ';', round(mean(auqct_het_rf/auqct_unif_rf)-1+1.96*sd(auqct_het_rf/auqct_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100,']'))

round(100*(mean(auqcs_het_rf/auqcs_unif_rf)-1),2)
print(paste0('[', round(mean(auqcs_het_rf/auqcs_unif_rf)-1-1.96*sd(auqcs_het_rf/auqcs_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100, ';', round(mean(auqcs_het_rf/auqcs_unif_rf)-1+1.96*sd(auqcs_het_rf/auqcs_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100,']'))

round(100*(mean(auqcx_het_rf/auqcx_unif_rf)-1),2)
print(paste0('[', round(mean(auqcx_het_rf/auqcx_unif_rf)-1-1.96*sd(auqcx_het_rf/auqcx_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100, ';', round(mean(auqcx_het_rf/auqcx_unif_rf)-1+1.96*sd(auqcx_het_rf/auqcx_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100,']'))

round(100*(mean(auqccf_het/auqccf_unif)-1),2)
print(paste0('[', round(mean(auqccf_het/auqccf_unif)-1-1.96*sd(auqccf_het/auqccf_unif)/sqrt(length(auqccf_unif)),4)*100, ';', round(mean(auqccf_het/auqccf_unif)-1+1.96*sd(auqccf_het/auqccf_unif)/sqrt(length(auqccf_unif)),4)*100,']'))

round(100*(mean(auqccb_het/auqccb_unif)-1),2)
print(paste0('[', round(mean(auqccb_het/auqccb_unif)-1-1.96*sd(auqccb_het/auqccb_unif)/sqrt(length(auqccb_unif)),4)*100, ';', round(mean(auqccb_het/auqccb_unif)-1+1.96*sd(auqccb_het/auqccb_unif)/sqrt(length(auqccb_unif)),4)*100,']'))

mean(auqct_het_rf)
mean(auqcs_het_rf)
mean(auqcx_het_rf)
mean(auqct_unif_rf)
mean(auqcs_unif_rf)
mean(auqcx_unif_rf)
mean(auqccb_unif)
mean(auqccb_het)
