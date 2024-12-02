library(MASS)
library(ranger)
library(data.table)
library(magrittr)
library(dplyr)
library(latex2exp)
library(grf)
library(bartCause)
library(isotone)

## required functions
source('code/helper/functions_simulation.R')
source('code/simulation_scenarios.R')
source('code/helper/functions_calibration.R')

## set seed
set.seed(12102024)

## get sampling parameters for each simulation scenario
for(i_scen in 1:3){
  
  rel <- get(paste0('rel', i_scen))
  treat <- get(paste0('treat', i_scen))
  scen <- paste0('sim_scen', i_scen)
  
  ## generate training data
  dt_train <- generate_data(p_n = 100000, p_rel = rel, p_treat = treat)
  dt_train_cal <- dt_train[1:80000,]
  dt_val <- dt_train[80001:90000,]
  dt_cal <- dt_train[-(1:90000),]
  
  ## train restricted outcome model
  rf_rest <- ranger(Y ~ ., data = dt_train[W==0, -c('W','mux','taux','expy'), with = F], max.depth = 5)
  
  ## tuning according to Brier Score
  max_depth_values <- c(5, 15, 25, 35, 45, 55)
  rf_tuned_list <- list()
  brier_score <- c()
  for(i_depth in 1:length(max_depth_values)){
    rf_tuned <- ranger(Y ~ ., data = dt_train_cal[W==0, -c('W','mux','taux','expy'), with = F], max.depth = max_depth_values[i_depth])
    rf_tuned_list[[i_depth]] <- rf_tuned
    brier_score <- c(brier_score, mean((dt_val[W==0,]$Y-predict(rf_tuned, dt_val[W==0,])$predictions)^2))
  }
  rf_best <- rf_tuned_list[[which.min(brier_score)]]
  
  
  # Platt calibrated outcome model
  Platt_outcome_func <- create_Platt_outcome_prediction(dt_cal[W == 0,], rf_best)
  
  ## isotonic regression outcome model 
  isotonic_outcome_func <- create_isotonic_outcome_prediction(dt_cal[W == 0,], rf_best)
  
  ############## Calculation of sampling parameters for all prediction methods ###############################
  dt_test <- generate_data(p_n = 100000, p_rel = rel, p_treat = treat)
  
  ## Results RF rest
  dt_test$pred_prob <- predict(rf_rest, data = dt_test)$predictions
  results_rf_rest <- get_sampling_parameters(copy(dt_test), rf_rest)
  
  ## Results RF Platt calibrated
  dt_test$pred_prob <- Platt_outcome_func(dt_test)
  results_rf_platt <- get_sampling_parameters(copy(dt_test), rf_best, p_outcome_func = Platt_outcome_func)
  
  ## Results RF isotonic regression calibrated
  dt_test$pred_prob <- isotonic_outcome_func(dt_test)
  results_rf_isotonic <- get_sampling_parameters(copy(dt_test), rf_best, p_outcome_func = isotonic_outcome_func)
  
  ## Results RF tuned
  dt_test$pred_prob <- predict(rf_best, data = dt_test)$predictions
  results_rf_tuned <- get_sampling_parameters(copy(dt_test), rf_best)
  
  ## save results
  saveRDS(results_rf_rest, paste0('results/sim_sampling_parameters/', scen, '/param_rf_rest.RDS'))
  saveRDS(results_rf_platt, paste0('results/sim_sampling_parameters/', scen, '/param_rf_platt.RDS'))
  saveRDS(results_rf_isotonic, paste0('results/sim_sampling_parameters/', scen, '/param_rf_isotonic.RDS'))
  saveRDS(results_rf_tuned, paste0('results/sim_sampling_parameters/', scen, '/param_rf_tuned.RDS'))
  
}

###### Variance quotient ######################################################
scenarios <- c('sim_scen1', 'sim_scen2', 'sim_scen3')

layout(matrix(c(1,2,3,4,4,4), ncol=3, byrow=TRUE), heights = c(0.9,0.1))
par(mar = c(3, 3, 3, 3))
for(i_scen in 1:length(scenarios)){
  results_rf_rest <- readRDS(paste0('results/sim_sampling_parameters/', scenarios[i_scen], '/param_rf_rest.RDS'))
  results_rf_platt <- readRDS(paste0('results/sim_sampling_parameters/', scenarios[i_scen], '/param_rf_platt.RDS'))
  results_rf_isotonic <- readRDS(paste0('results/sim_sampling_parameters/', scenarios[i_scen], '/param_rf_isotonic.RDS'))
  results_rf_tuned <- readRDS(paste0('results/sim_sampling_parameters/', scenarios[i_scen], '/param_rf_tuned.RDS'))
  
  
  plot(results_rf_rest$pH, (results_rf_rest$Q_V_est-results_rf_rest$Q_V_ass)/results_rf_rest$Q_V_ass*100, type = 'l', xlab = TeX('$p_H$'), 
       ylab = TeX('% estimation error for $Q_V$'), main = paste0('Simulation scenario ', i_scen), lwd = 3, ylim = c(-200, 400), mgp=c(2,1,0))
  lines(results_rf_rest$pH, (results_rf_platt$Q_V_est-results_rf_platt$Q_V_ass)/results_rf_platt$Q_V_ass*100, lty = 1, lwd = 3, col = 'blue')
  lines(results_rf_rest$pH, (results_rf_isotonic$Q_V_est-results_rf_isotonic$Q_V_ass)/results_rf_isotonic$Q_V_ass*100, lty = 1, col = 'green', lwd = 3)
  lines(results_rf_rest$pH, (results_rf_tuned$Q_V_est-results_rf_tuned$Q_V_ass)/results_rf_tuned$Q_V_ass*100, lty = 1, col = 'orange', lwd = 3)
  lines(results_rf_rest$pH, rep(0, length(results_rf_rest$pH)), col = 'red', lty = 3)
}

##legend
par(mar=c(0, 0, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('red','black', 'blue', 'green', 'orange')
legend(x = "bottom",inset = 0,
       legend = c(TeX('Optimal'), 'HS rest', 'HS Platt', 'HS isotonic', 'HS tuned'), 
       col=plot_colors, lty = c(3,1,1,1,1), lwd=c(3,3), cex=1, horiz = TRUE, bty = "n")