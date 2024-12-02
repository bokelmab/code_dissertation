library(MASS)
library(ranger)
library(data.table)
library(magrittr)
library(dplyr)
#library(latex2exp)
#library(grf)
#library(bartCause)
#library(isotone)

## required functions
source('code/helper/functions_simulation.R')
source('code/simulation_scenarios.R')

## set seed
set.seed(12102024)

## get sampling parameters for each simulation scenario
for(i_scen in 1:3){
  
  rel <- get(paste0('rel', i_scen))
  treat <- get(paste0('treat', i_scen))
  scen <- paste0('sim_scen', i_scen)
  
  ## generate training data
  dt_train <- generate_data(p_n = 100000, p_rel = rel, p_treat = treat)
  dt_train_small <- copy(dt_train[1:round(nrow(dt_train)*0.1),])
  
  ## train restricted outcome model
  rf_small <- ranger(Y ~ ., data = dt_train_small[W==0, -c('W','mux','taux','expy'), with = F], max.depth = 5)
  
  ############## Calculation of sampling parameters for all prediction methods ###############################
  dt_test <- generate_data(p_n = 100000, p_rel = rel, p_treat = treat)
  
  ## Results RF rest
  dt_test$pred_prob <- predict(rf_small, data = dt_test)$predictions
  results_rf_small <- get_sampling_parameters(copy(dt_test), rf_small)
  
  ## save results
  saveRDS(results_rf_small, paste0('results/sim_sampling_parameters/', scen, '/param_rf_small.RDS'))
  
}

