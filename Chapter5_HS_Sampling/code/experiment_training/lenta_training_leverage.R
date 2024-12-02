library(MASS)
library(ranger)
library(data.table)
library(magrittr)
library(latex2exp)
library(grf)
#library(bartCause)

## required functions
source('code/helper/functions_simulation.R')
source('code/simulation_scenarios.R')
source('code/helper/functions_effect_estimation.R')

## choice of simulation setting and sampling parameters
samp_param_method <- 'leverage'

## read and prepare Lenta data
dt_exp <- fread(paste0('results/sim_sampling_parameters/', 'lenta', '/dt_exp.csv'))

prob_W <- 0.75

y_data <- dt_exp$Y
dt_exp <- as.data.table(model.matrix(Y~-1+., dt_exp)) ## for causal random forest
dt_exp$Y <- y_data

auqct_unif_glm <- c()
auqct_het_glm <- c()
auqct_unif_rf <- c()
auqcs_unif_rf <- c()
auqcx_unif_rf <- c()
auqct_het_rf <- c()
auqcs_het_rf <- c()
auqcx_het_rf <- c()
auqccf_het <- c()
auqccf_unif <- c()
for(i_sim in 1:1000){
  
  idx_train <- sample(1:nrow(dt_exp), 100000)
  dt_train <- dt_exp[idx_train,]
  idx_unif <- sample(1:nrow(dt_train), 20000)
  ## leverage score sampling probabilites
  X_matr <- model.matrix(Y~-1+., dt_train[,-c('W'), with = F])
  H <- X_matr %*% solve(t(X_matr) %*% X_matr) %*% t(X_matr)
  lev_scores <- diag(H)
  idx_lev <- sample(1:nrow(dt_train), 20000, prob = lev_scores/sum(lev_scores))
  train_unif <- dt_train[idx_unif,]
  train_het <- dt_train[idx_lev,]
  dt_test <- dt_exp[-idx_train,]
  
  ############ Training on completely randomly sampled data ##############################
  
  ## T-learner GLM
  rf0_unif <- glm(Y ~ ., data = train_unif[W==0,-c('W', 'var_group'), with = F], family = 'binomial')
  rf1_unif <- glm(Y ~ ., data = train_unif[W==1,-c('W', 'var_group'), with = F], family = 'binomial')
  predt_unif <- predict(rf1_unif, newdata = dt_test, type = 'response')-predict(rf0_unif, newdata = dt_test, type = 'response')
  auqct_unif_glm <- c(auqct_unif_glm, calc_auq(predt_unif, dt_test))
  
  ## T-learner RF
  rf0_unif <- ranger(Y ~ ., data = train_unif[W==0,-c('W', 'var_group'), with = F])
  rf1_unif <- ranger(Y ~ ., data = train_unif[W==1,-c('W', 'var_group'), with = F])
  predt_unif <- predict(rf1_unif, data = dt_test)$predictions-predict(rf0_unif, data = dt_test)$predictions
  auqct_unif_rf <- c(auqct_unif_rf, calc_auq(predt_unif, dt_test))
  
  ## S-Learner
  rfs_unif <- ranger(Y ~ ., data = train_unif[,-c('var_group'), with = F])
  preds_unif <- predict(rfs_unif, data = copy(dt_test)[,W:=1])$predictions-predict(rfs_unif, data = copy(dt_test)[,W:=0])$predictions
  auqcs_unif_rf <- c(auqcs_unif_rf, calc_auq(preds_unif, dt_test))
  
  ## X-Learner
  train_unif_x <- copy(train_unif[,-c('var_group'), with = F])
  train_unif_x1 <- train_unif_x[W==1,]
  train_unif_x0 <- train_unif_x[W==0,]
  train_unif_x1$Y <- train_unif_x1$Y - predict(rf0_unif, train_unif_x1[,-c('W'),with = F])$predictions
  train_unif_x0$Y <- predict(rf1_unif, train_unif_x0[,-c('W'),with = F])$predictions - train_unif_x0$Y
  rfX0_unif <- ranger(Y ~ ., data = train_unif_x0[,-c('W'), with = F])
  rfX1_unif <- ranger(Y ~ ., data = train_unif_x1[,-c('W'), with = F])
  predx_unif <- 0.25*predict(rfX0_unif, data = dt_test)$predictions+0.75*predict(rfX1_unif, data = dt_test)$predictions
  auqcx_unif_rf <- c(auqcx_unif_rf, calc_auq(predx_unif, dt_test))
  
  ## Causal random forest
  cf_unif <- causal_forest(X = train_unif[,-c('W', 'Y','var_group'), with = F], Y = train_unif$Y, W = train_unif$W, W.hat = 0.75, num.trees = 500)
  predcf_unif <- predict(cf_unif, dt_test[,-c('W', 'Y','var_group'), with = F])$predictions
  ecdfcf_unif <- ecdf(predcf_unif)
  auqccf_unif <- c(auqccf_unif, calc_auq(predcf_unif, dt_test))
  
  ## causal bart
  # cbart <- bartc(response = train_unif$Y, treatment = train_unif$W, 
  #                confounders = train_unif[,-c('W', 'Y','var_group'), with = F], 
  #                n.samples = 100L, n.burn = 15L, n.chains = 2L, method.trt = rep(prob_W, nrow(train_unif)),keepTrees = T)
  # predccb_unif <- apply(predict(cbart, newdata = copy(dt_test[, ps:= 1]), type = 'icate'),2,mean)
  # ecdfcb_unif <- ecdf(predccb_unif)
  # auqccb_unif <- c(auqccb_unif, calc_auq(predccb_unif, dt_test))
  
  ########## Training on leverage sampled data ##################
  
  ## T-learner GLM
  rf0_het <- glm(Y ~ ., data = train_het[W==0,-c('W', 'var_group'), with = F], family = 'binomial') 
  rf1_het <- glm(Y ~ ., data = train_het[W==1,-c('W', 'var_group'), with = F], family = 'binomial')
  predt_het <- predict(rf1_het, newdata = dt_test, type = 'response')-predict(rf0_het, newdata = dt_test, type = 'response')
  auqct_het_glm <- c(auqct_het_glm, calc_auq(predt_het, dt_test))
  
  ## T-learner
  rf0_het <- ranger(Y ~ ., data = train_het[W==0,-c('W', 'var_group'), with = F])
  rf1_het <- ranger(Y ~ ., data = train_het[W==1,-c('W', 'var_group'), with = F])
  predt_het <- predict(rf1_het, data = dt_test)$predictions-predict(rf0_het, data = dt_test)$predictions
  auqct_het_rf <- c(auqct_het_rf, calc_auq(predt_het, dt_test))
  
  ## S-Learner
  rfs_het <- ranger(Y ~ ., data = train_het[,-c('var_group'), with = F])
  preds_het <- predict(rfs_het, data = copy(dt_test)[,W:=1])$predictions-predict(rfs_het, data = copy(dt_test)[,W:=0])$predictions
  auqcs_het_rf <- c(auqcs_het_rf, calc_auq(preds_het, dt_test))
  
  ## X-Learner
  train_het_x <- copy(train_het[,-c('var_group'), with = F])
  train_het_x1 <- train_het_x[W==1,]
  train_het_x0 <- train_het_x[W==0,]
  train_het_x1$Y <- train_het_x1$Y - predict(rf0_het, train_het_x1[,-c('W'),with = F])$predictions
  train_het_x0$Y <- predict(rf1_het, train_het_x0[,-c('W'),with = F])$predictions - train_het_x0$Y
  rfX0_het <- ranger(Y ~ ., data = train_het_x0[,-c('W'), with = F])
  rfX1_het <- ranger(Y ~ ., data = train_het_x1[,-c('W'), with = F])
  predx_het <- 0.25*predict(rfX0_het, data = dt_test)$predictions+0.75*predict(rfX1_het, data = dt_test)$predictions
  auqcx_het_rf <- c(auqcx_het_rf, calc_auq(predx_het, dt_test))
  
  ## Causal random forest
  cf_het <- causal_forest(X = train_het[,-c('W', 'Y','var_group'), with = F], Y = train_het$Y, W = train_het$W, W.hat = 0.75, num.trees = 500)
  predcf_het <- predict(cf_het, dt_test[,-c('W', 'Y','var_group'), with = F])$predictions
  ecdfcf_het <- ecdf(predcf_het)
  auqccf_het <- c(auqccf_het, calc_auq(predcf_het, dt_test))
  
  ## Causal Bart
  # cbart_het <- bartc(response = train_het$Y, treatment = train_het$W, 
  #                    confounders = train_het[,-c('W', 'Y','var_group'), with = F], 
  #                    n.samples = 100L, n.burn = 15L, n.chains = 2L, method.trt = rep(prob_W, nrow(train_het)),keepTrees = T)
  # predccb_het <- apply(predict(cbart_het, newdata = copy(dt_test[, ps:= 1]), type = 'icate'),2,mean)
  # ecdfcb_het <- ecdf(predccb_het)
  # auqccb_het <- c(auqccb_het, calc_auq(predccb_het, dt_test))
  
  print(i_sim)
  ## evaluation T-learner glm
  100*(mean(auqct_het_glm/auqct_unif_glm)-1)
  print(paste0('[', round(mean(auqct_het_glm/auqct_unif_glm)-1-1.96*sd(auqct_het_glm/auqct_unif_glm)/sqrt(length(auqct_unif_glm)),4)*100, ';', round(mean(auqct_het_glm/auqct_unif_glm)-1+1.96*sd(auqct_het_glm/auqct_unif_glm)/sqrt(length(auqct_unif_glm)),4)*100,']'))
  
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
  # 100*(mean(auqccb_het/auqccb_unif)-1)
  # print(paste0('[', round(mean(auqccb_het/auqccb_unif)-1-1.96*sd(auqccb_het/auqccb_unif)/sqrt(length(auqccb_het)),4)*100, ';', round(mean(auqccb_het/auqccb_unif)-1+1.96*sd(auqccb_het/auqccb_unif)/sqrt(length(auqccb_het)),4)*100,']'))
  # 
  
  print(mean(auqct_het_glm))
  print(mean(auqct_het_rf))
  print(mean(auqcs_het_rf))
  print(mean(auqcx_het_rf))
  print(mean(auqccf_het))
  # print(mean(auqccb_het))
  
  ## save results
  saveRDS(auqct_het_glm, paste0('results/', 'lenta/', samp_param_method, '/auqct_het_glm.RDS'))
  saveRDS(auqct_het_rf, paste0('results/', 'lenta/', samp_param_method, '/auqct_het_rf.RDS'))
  saveRDS(auqcs_het_rf, paste0('results/', 'lenta/', samp_param_method, '/auqcs_het_rf.RDS'))
  saveRDS(auqcx_het_rf, paste0('results/', 'lenta/', samp_param_method, '/auqcx_het_rf.RDS'))
  saveRDS(auqccf_het, paste0('results/', 'lenta/', samp_param_method, '/auqccf_het.RDS'))
  # saveRDS(auqccb_het, paste0('results/', 'lenta', '/auqccb_het.RDS'))
  
  
  saveRDS(auqct_unif_glm, paste0('results/', 'lenta/', samp_param_method, '/auqct_unif_glm.RDS'))
  saveRDS(auqct_unif_rf, paste0('results/', 'lenta/', samp_param_method, '/auqct_unif_rf.RDS'))
  saveRDS(auqcs_unif_rf, paste0('results/', 'lenta/', samp_param_method, '/auqcs_unif_rf.RDS'))
  saveRDS(auqcx_unif_rf, paste0('results/', 'lenta/', samp_param_method, '/auqcx_unif_rf.RDS'))
  saveRDS(auqccf_unif, paste0('results/', 'lenta/', samp_param_method, '/auqccf_unif.RDS'))
  # saveRDS(auqccb_unif, paste0('results/', 'lenta', '/auqccb_unif.RDS'))
  
}

## calculate results for paper
round(100*(mean(auqccf_het/auqccf_unif)-1),2)
print(paste0('[', round(mean(auqccf_het/auqccf_unif)-1-1.96*sd(auqccf_het/auqccf_unif)/sqrt(length(auqccf_unif)),4)*100, ';', round(mean(auqccf_het/auqccf_unif)-1+1.96*sd(auqccf_het/auqccf_unif)/sqrt(length(auqccf_unif)),4)*100,']'))

round(100*(mean(auqct_het_glm/auqct_unif_glm)-1),2)
print(paste0('[', round(mean(auqct_het_glm/auqct_unif_glm)-1-1.96*sd(auqct_het_glm/auqct_unif_glm)/sqrt(length(auqct_unif_glm)),4)*100, ';', round(mean(auqct_het_glm/auqct_unif_glm)-1+1.96*sd(auqct_het_glm/auqct_unif_glm)/sqrt(length(auqct_unif_glm)),4)*100,']'))

round(100*(mean(auqct_het_rf/auqct_unif_rf)-1),2)
print(paste0('[', round(mean(auqct_het_rf/auqct_unif_rf)-1-1.96*sd(auqct_het_rf/auqct_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100, ';', round(mean(auqct_het_rf/auqct_unif_rf)-1+1.96*sd(auqct_het_rf/auqct_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100,']'))

round(100*(mean(auqcs_het_rf/auqcs_unif_rf)-1),2)
print(paste0('[', round(mean(auqcs_het_rf/auqcs_unif_rf)-1-1.96*sd(auqcs_het_rf/auqcs_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100, ';', round(mean(auqcs_het_rf/auqcs_unif_rf)-1+1.96*sd(auqcs_het_rf/auqcs_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100,']'))

round(100*(mean(auqcx_het_rf/auqcx_unif_rf)-1),2)
print(paste0('[', round(mean(auqcx_het_rf/auqcx_unif_rf)-1-1.96*sd(auqcx_het_rf/auqcx_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100, ';', round(mean(auqcx_het_rf/auqcx_unif_rf)-1+1.96*sd(auqcx_het_rf/auqcx_unif_rf)/sqrt(length(auqct_unif_rf)),4)*100,']'))


# round(100*(mean(auqccb_het/auqccb_unif)-1),2)
# print(paste0('[', round(mean(auqccb_het/auqccb_unif)-1-1.96*sd(auqccb_het/auqccb_unif)/sqrt(length(auqccb_unif)),4)*100, ';', round(mean(auqccb_het/auqccb_unif)-1+1.96*sd(auqccb_het/auqccb_unif)/sqrt(length(auqccb_unif)),4)*100,']'))


mean(auqct_het_glm)
mean(auqct_unif_glm)
mean(auqct_het_rf)
mean(auqcs_het_rf)
mean(auqcx_het_rf)
mean(auqct_unif_rf)
mean(auqcs_unif_rf)
mean(auqcx_unif_rf)
mean(auqccf_unif)
mean(auqccf_het)