library(MASS)
library(ranger)
library(data.table)
library(magrittr)
#library(latex2exp)
#library(grf)
#library(bartCause)

## required functions
source('code/helper/functions_simulation.R')
source('code/simulation_scenarios.R')


for(i_scen in 1:3){
  
  scen_choice <- i_scen
  
  ## evaluate each method to obtain sampling parameters
  samp_methods <- c('small')
  for(i_samp_met in 1:length(samp_methods)){
    
    samp_param_method <- samp_methods[i_samp_met]
    print(samp_param_method)  
    
    ## create results folder
    path_results <- paste0('results/sim_scen', scen_choice, '/', samp_param_method)
    dir.create(path_results, path_results)
    
    ## simulation scenario
    rel <- get(paste0('rel', scen_choice))
    treat <- get(paste0('treat', scen_choice))
    scen <- paste0('sim_scen', scen_choice)
    
    ## read sampling parameters
    samp_param <- readRDS(paste0('results/sim_sampling_parameters/', scen, '/param_rf_', samp_param_method, '.RDS'))
    outcome_model <- samp_param$outcome_model
    pH <- samp_param$pH_choose
    pL <- 1-pH
    RH <- samp_param$R_H_choose
    share_high <- RH*(1-pL)
    saveRDS(list(pH = pH, RH = RH), paste0(path_results, '/sampling_parameters.RDS'))
    
    ########## Evaluation part #####################################################
    data_train <- generate_data(p_n = 50000, p_rel = rel, p_treat = treat)
    data_complete <- generate_data(p_n = 250000, p_rel = rel, p_treat = treat)
    
    rf0_unif <- ranger(Y ~ ., data = data_train[W==0,-c('W','mux','taux','expy'), with = F])
    rf1_unif <- ranger(Y ~ ., data = data_train[W==1,-c('W','mux','taux','expy'), with = F]) 
    data_complete$predt <- predict(rf1_unif, data = data_complete)$predictions-predict(rf0_unif, data = data_complete)$predictions
    data_complete$pred_prob <- samp_param$outcome_func(data_complete) 
    
    data_complete$var_group <- ifelse(data_complete$pred_prob > quantile(data_complete$pred_prob,pL), 2, 1)
    fwrite(data_complete, paste0(path_results, '/data_complete.csv'))
    
    est_ate_unif <- c()
    est_ate_het <- c()
    est_ate_ca <- c()
    est_ate_hsca <- c()
    qini_table_het <- NULL
    qini_table_unif <- NULL
    qini_table_ca <- NULL
    qini_table_hsca <- NULL
    
    for(i_sim in 1:1000){
      
      data <- generate_data(p_n = 700000, p_rel = rel, p_treat = treat)
      data$pred_prob <- samp_param$outcome_func(data) 
      data$var_group <- ifelse(data$pred_prob > quantile(data_complete$pred_prob,pL), 2, 1)
      
      test_unif <- data[sample(1:nrow(data), 10000),]
      idx_var_group1 <- sample(which(data$var_group == 1), round((1-share_high)*10000))
      idx_var_group2 <- sample(which(data$var_group == 2), round(share_high*10000))
      test_het <- data[c(idx_var_group1,idx_var_group2),]
      
      test_unif$predt <- predict(rf1_unif, data = test_unif)$predictions-predict(rf0_unif, data = test_unif)$predictions
      test_het$predt <- predict(rf1_unif, data = test_het)$predictions-predict(rf0_unif, data = test_het)$predictions
      
      ## ATE treatment prio
      qini_unif <- c()
      qini_het <- c()
      qini_ca <- c()
      qini_hsca <- c()
      for(i_cut in 9:0){
        
        cutoff <- i_cut/10
        quant <- quantile(data_complete$predt, cutoff)
        
        ## ate estimation for t-learner on uniform data
        ate_unif <- mean(test_unif[predt >= quant & W == 1,]$Y)-mean(test_unif[predt >= quant & W == 0,]$Y)
        ate_ca <- mean(test_unif[predt >= quant & W == 1,]$Y-test_unif[predt >= quant & W == 1,]$mux-0.5*test_unif[predt >= quant & W == 1,]$taux)-mean(test_unif[predt >= quant & W == 0,]$Y-test_unif[predt >= quant & W == 0,]$mux-0.5*test_unif[predt >= quant & W == 0,]$taux)
        weight_group1 <- mean(data_complete[predt >= quant,]$var_group == 1)
        ate_v1 <- mean(test_het[predt >= quant & W == 1 & var_group == 1,]$Y)-mean(test_het[predt >= quant & W == 0 & var_group == 1,]$Y)
        ate_v2 <- mean(test_het[predt >= quant & W == 1 & var_group == 2,]$Y)-mean(test_het[predt >= quant & W == 0 & var_group == 2,]$Y)
        ate_het <- weight_group1*ate_v1 + (1-weight_group1)*ate_v2
        ate_v1_ca <- mean(test_het[predt >= quant & W == 1 & var_group == 1,]$Y-test_het[predt >= quant & W == 1 & var_group == 1,]$mux-0.5*test_het[predt >= quant & W == 1 & var_group == 1,]$taux)-mean(test_het[predt >= quant & W == 0 & var_group == 1,]$Y-test_het[predt >= quant & W == 0 & var_group == 1,]$mux-0.5*test_het[predt >= quant & W == 0 & var_group == 1,]$taux)
        ate_v2_ca <- mean(test_het[predt >= quant & W == 1 & var_group == 2,]$Y-test_het[predt >= quant & W == 1 & var_group == 2,]$mux-0.5*test_het[predt >= quant & W == 1 & var_group == 2,]$taux)-mean(test_het[predt >= quant & W == 0 & var_group == 2,]$Y-test_het[predt >= quant & W == 0 & var_group == 2,]$mux-0.5*test_het[predt >= quant & W == 0 & var_group == 2,]$taux)
        ate_hsca <- weight_group1*ate_v1_ca + (1-weight_group1)*ate_v2_ca
        rm(weight_group1, ate_v1, ate_v2)
        
        qini_unif <- c(qini_unif, ate_unif*(10-i_cut))
        qini_het <- c(qini_het, ate_het*(10-i_cut))
        qini_ca <- c(qini_ca, ate_ca*(10-i_cut))
        qini_hsca <- c(qini_hsca, ate_hsca*(10-i_cut))
        
      }
      
      qini_table_unif %<>% rbind(qini_unif)
      qini_table_het %<>% rbind(qini_het)
      qini_table_ca %<>% rbind(qini_ca)
      qini_table_hsca %<>% rbind(qini_hsca)
      
      
      ## ate est
      est_ate_unif <- c(est_ate_unif, mean(test_unif[W==1,]$Y)-mean(test_unif[W==0,]$Y))
      ate1 <- mean(test_het[W==1 & var_group == 1,]$Y)-mean(test_het[W==0 & var_group == 1,]$Y)
      ate2 <- mean(test_het[W==1 & var_group == 2,]$Y)-mean(test_het[W==0 & var_group == 2,]$Y)
      est_ate_het <- c(est_ate_het, pL*ate1+(1-pL)*ate2)
      est_ate_ca <- c(est_ate_ca, mean(test_unif[W==1,]$Y-test_unif[W==1,]$mux-0.5*test_unif[W==1,]$taux)-mean(test_unif[W==0,]$Y-test_unif[W==0,]$mux-0.5*test_unif[W==0,]$taux))
      ate1_ca <- mean(test_het[W==1 & var_group == 1,]$Y-test_het[W==1 & var_group == 1,]$mux-0.5*test_het[W==1 & var_group == 1,]$taux)-mean(test_het[W==0 & var_group == 1,]$Y-test_het[W==0 & var_group == 1,]$mux-0.5*test_het[W==0 & var_group == 1,]$taux)
      ate2_ca <- mean(test_het[W==1 & var_group == 2,]$Y-test_het[W==1 & var_group == 2,]$mux-0.5*test_het[W==1 & var_group == 2,]$taux)-mean(test_het[W==0 & var_group == 2,]$Y-test_het[W==0 & var_group == 2,]$mux-0.5*test_het[W==0 & var_group == 2,]$taux)
      est_ate_hsca <- c(est_ate_hsca, pL*ate1_ca+(1-pL)*ate2_ca)
      print(i_sim)
      
      if(i_sim %% 100 == 0){
        print(1-var(est_ate_het)/var(est_ate_unif))
        print(var(est_ate_het))
        print(samp_param$calc_true_var_ate(samp_param$pH_choose, samp_param$R_H_choose))
        print(var(est_ate_unif))
        print(samp_param$var_rand)
      }
      
      plot(0:10, c(0,apply(qini_table_unif,2,mean)),type = 'l', ylim = c(0, 12*mean(data_complete$taux)))
      lines(0:10, c(0,apply(qini_table_het,2,mean)), col = 'blue')
      lines(0:10, c(0,apply(qini_table_ca,2,mean)), col = 'red')
      lines(0:10, c(0,apply(qini_table_hsca,2,mean)), col = 'green')
      
      
    }
    round((1-var(est_ate_ca)/var(est_ate_unif))*100,1)
    round((1-var(est_ate_het)/var(est_ate_unif))*100,1)
    round((1-var(est_ate_hsca)/var(est_ate_unif))*100,1)
    
    round((1-apply(qini_table_ca,2,var)/apply(qini_table_unif,2,var))*100,1)
    round((1-apply(qini_table_het,2,var)/apply(qini_table_unif,2,var))*100,1)
    round((1-apply(qini_table_hsca,2,var)/apply(qini_table_unif,2,var))*100,1)
    
    ## save results
    saveRDS(est_ate_ca, paste0(path_results, '/est_ate_ca.RDS'))
    saveRDS(est_ate_het, paste0(path_results, '/est_ate_het.RDS'))
    saveRDS(est_ate_hsca, paste0(path_results, '/est_ate_hsca.RDS'))
    saveRDS(est_ate_unif, paste0(path_results, '/est_ate_unif.RDS'))
    
    saveRDS(qini_table_ca, paste0(path_results, '/qini_table_ca.RDS'))
    saveRDS(qini_table_het, paste0(path_results, '/qini_table_het.RDS'))
    saveRDS(qini_table_hsca, paste0(path_results, '/qini_table_hsca.RDS'))
    saveRDS(qini_table_unif, paste0(path_results, '/qini_table_unif.RDS'))
  }
}
