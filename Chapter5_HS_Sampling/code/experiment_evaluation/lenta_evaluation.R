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


## read experiment data
dt_exp <- fread(paste0('results/sim_sampling_parameters/', 'lenta', '/dt_exp.csv'))
prob_W <- 0.75

## evaluate each method to obtain sampling parameters
samp_methods <- c('isotonic', 'platt', 'rest', 'tuned', 'small')
for(i_samp_met in 1:length(samp_methods)){
  
  samp_param_method <- samp_methods[i_samp_met]
  print(samp_param_method)  
  
  ## create results folder
  path_results <- paste0('results/lenta/', samp_param_method)
  dir.create(path_results)
  
  ## read sampling parameters
  samp_param <- readRDS(paste0('results/sim_sampling_parameters/lenta/param_rf_', samp_param_method, '.RDS'))
  outcome_model <- samp_param$outcome_model
  pH <- samp_param$pH_choose
  pL <- 1-pH
  RH <- samp_param$R_H_choose
  share_high <- RH*(1-pL)
  saveRDS(list(pH = pH, RH = RH), paste0(path_results, '/sampling_parameters.RDS'))
  
  ########## Evaluation part #####################################################
  ## define strata in experimental data
  pred_prob <- samp_param$outcome_func(dt_exp)
  dt_exp$var_group <- ifelse(pred_prob > quantile(pred_prob,pL), 2, 1)
  
  ## training of uplift model
  idx_train <- sample(1:nrow(dt_exp), 100000, replace = F)
  dt_train <- dt_exp[idx_train,]
  dt_test <- dt_exp[-idx_train,]
  idx_val <- sample(1:100000, 10000)
  
  ## tuning according to Brier Score for covariate adjustment
  max_depth_values <- c(5, 15, 25, 35, 45, 55)
  brier_score <- c()
  for(i_depth in 1:length(max_depth_values)){
    rf_tuned <- ranger(Y ~ ., data = dt_train[-idx_val,][W==0, -c('W'), with = F], max.depth = max_depth_values[i_depth])
    brier_score <- c(brier_score, mean((dt_train[idx_val,][W==0,]$Y-predict(rf_tuned, dt_train[idx_val,][W==0,])$predictions)^2))
  }
  rf0 <- ranger(Y ~ ., data = dt_train[W==0,-c('W','var_group'), with = F], max.depth = max_depth_values[which.min(brier_score)])
  rf1 <- ranger(Y ~ ., data = dt_train[W==1,-c('W','var_group'), with = F], max.depth = max_depth_values[which.min(brier_score)])
  
  ## make predictions on test set
  dt_test$predt <- predict(rf1, dt_test)$predictions - predict(rf0, dt_test)$predictions
  dt_test$phi_hat <- prob_W*predict(rf0, dt_test)$predictions + (1-prob_W)*predict(rf1, dt_test)$predictions
  
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
    idx_var_group2 <- sample(which(dt_test$var_group == 2), 10000*share_high)
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
    
    ## plot Qini curve
    plot(0:10, c(0,apply(qini_table_unif,2,mean)),type = 'l')
    lines(0:10, c(0,apply(qini_table_het,2,mean)), col = 'blue')
    lines(0:10, c(0,apply(qini_table_ca,2,mean)), col = 'red')
    lines(0:10, c(0,apply(qini_table_hsca,2,mean)), col = 'green')
    
    ## save results
    saveRDS(est_ate_ca, paste0('results/', 'lenta/', samp_param_method, '/est_ate_ca.RDS'))
    saveRDS(est_ate_het, paste0('results/', 'lenta/', samp_param_method, '/est_ate_het.RDS'))
    saveRDS(est_ate_hsca, paste0('results/', 'lenta/', samp_param_method, '/est_ate_hsca.RDS'))
    saveRDS(est_ate_unif, paste0('results/', 'lenta/', samp_param_method, '/est_ate_unif.RDS'))
    
    saveRDS(var_ate_ca, paste0('results/', 'lenta/', samp_param_method, '/var_ate_ca.RDS'))
    saveRDS(var_ate_het, paste0('results/', 'lenta/', samp_param_method, '/var_ate_het.RDS'))
    saveRDS(var_ate_hsca, paste0('results/', 'lenta/', samp_param_method, '/var_ate_hsca.RDS'))
    saveRDS(var_ate_unif, paste0('results/', 'lenta/', samp_param_method, '/var_ate_unif.RDS'))
    
    saveRDS(qini_table_ca, paste0('results/', 'lenta/', samp_param_method, '/qini_table_ca.RDS'))
    saveRDS(qini_table_het, paste0('results/', 'lenta/', samp_param_method, '/qini_table_het.RDS'))
    saveRDS(qini_table_hsca, paste0('results/', 'lenta/', samp_param_method, '/qini_table_hsca.RDS'))
    saveRDS(qini_table_unif, paste0('results/', 'lenta/', samp_param_method, '/qini_table_unif.RDS'))
    
    saveRDS(var_qini_table_ca, paste0('results/', 'lenta/', samp_param_method, '/var_qini_table_ca.RDS'))
    saveRDS(var_qini_table_het, paste0('results/', 'lenta/', samp_param_method, '/var_qini_table_het.RDS'))
    saveRDS(var_qini_table_hsca, paste0('results/', 'lenta/', samp_param_method, '/var_qini_table_hsca.RDS'))
    saveRDS(var_qini_table_unif, paste0('results/', 'lenta/', samp_param_method, '/var_qini_table_unif.RDS'))
    
    print(round((1-mean(var_ate_ca/var_ate_unif))*100,2))
    print(round((1-mean(var_ate_het/var_ate_unif))*100,2))
    print(round((1-mean(var_ate_hsca/var_ate_unif))*100,2))
  }
}
  
  plot(0:10, c(0, apply(qini_table_unif,2,mean)), type = 'l')
  lines(0:10, c(0, apply(qini_table_het,2,mean)), col = 'red')
  lines(0:10, c(0, apply(qini_table_ca,2,mean)), col = 'blue')
  lines(0:10, c(0, apply(qini_table_hsca,2,mean)), col = 'green')
  
  round((1-apply(var_qini_table_ca,2,mean)/apply(var_qini_table_unif,2,mean))*100,1)
  round((1-apply(var_qini_table_het,2,mean)/apply(var_qini_table_unif,2,mean))*100,1)
  round((1-apply(var_qini_table_hsca,2,mean)/apply(var_qini_table_unif,2,mean))*100,1)
  
  round((1-mean(var_ate_ca)/mean(var_ate_unif))*100,1)
  round((1-mean(var_ate_het)/mean(var_ate_unif))*100,1)
  round((1-mean(var_ate_hsca)/mean(var_ate_unif))*100,1)
