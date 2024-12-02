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
source('code/helper/functions_calibration.R')

## set seed
set.seed(12102024)

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
dt_train <- data[idx_base,]
dt_train_small <- copy(dt_train[1:round(nrow(dt_train)*0.1),])
dt_train_cal <- dt_train[1:round(nrow(dt_train)*0.8),]
dt_val <- dt_train[(round(nrow(dt_train)*0.8)+1):round(nrow(dt_train)*0.9),]
dt_cal <- dt_train[-(1:round(nrow(dt_train)*0.9)),]
dt_exp <- data[-idx_base,]
dt_test <- copy(dt_exp[sample(1:nrow(dt_exp), 100000),])

## train restricted outcome model
rf_rest <- ranger(Y ~ ., data = dt_train[W==0, -c('W'), with = F], max.depth = 5)

## train model with small data set
rf_small <- ranger(Y ~ ., data = dt_train_small[W==0, -c('W'), with = F], max.depth = 5)

## tuning according to Brier Score
max_depth_values <- c(5, 15, 25, 35, 45, 55)
rf_tuned_list <- list()
brier_score <- c()
for(i_depth in 1:length(max_depth_values)){
  rf_tuned <- ranger(Y ~ ., data = dt_train_cal[W==0, -c('W'), with = F], max.depth = max_depth_values[i_depth])
  rf_tuned_list[[i_depth]] <- rf_tuned
  brier_score <- c(brier_score, mean((dt_val[W==0,]$Y-predict(rf_tuned, dt_val[W==0,])$predictions)^2))
}
rf_best <- rf_tuned_list[[which.min(brier_score)]]

# define Platt calibrated outcome model with closure
Platt_outcome_func <- create_Platt_outcome_prediction(dt_cal[W == 0,], rf_best)

## define isotonic regression outcome model with closure
isotonic_outcome_func <- create_isotonic_outcome_prediction(dt_cal[W == 0,], rf_best)

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

## Results RF small
dt_test$pred_prob <- predict(rf_small, data = dt_test)$predictions
results_rf_small <- get_sampling_parameters(copy(dt_test), rf_small)

## save results
saveRDS(results_rf_rest, paste0('results/sim_sampling_parameters/', 'criteo', '/param_rf_rest.RDS'))
saveRDS(results_rf_platt, paste0('results/sim_sampling_parameters/', 'criteo', '/param_rf_platt.RDS'))
saveRDS(results_rf_isotonic, paste0('results/sim_sampling_parameters/', 'criteo', '/param_rf_isotonic.RDS'))
saveRDS(results_rf_tuned, paste0('results/sim_sampling_parameters/', 'criteo', '/param_rf_tuned.RDS'))
saveRDS(results_rf_small, paste0('results/sim_sampling_parameters/', 'criteo', '/param_rf_small.RDS'))
fwrite(dt_exp, paste0('results/sim_sampling_parameters/', 'criteo', '/dt_exp.csv'))
