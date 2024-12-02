##required libraries
library(data.table)
library(dplyr)
library(magrittr)
library(ranger)
library(grf)
#library(bartCause)
library(latex2exp)

## required source files
source('code/helper/functions_simulation.R')
source('code/helper/functions_calibration.R')

## set seed
set.seed(12102024)

## read and prepare Lenta data
data <- fread('data/lenta_dataset.csv')
data$W <- ifelse(data$group == 'test', 1, 0)
data$Y <- data$response_att
data$group <- NULL
data$response_att <- NULL

## impute important features
data[which(is.na(response_sms)), response_sms := mean(data$response_sms, na.rm = T)]
data[which(is.na(response_viber)), response_viber := mean(data$response_viber, na.rm = T)]
data[which(is.na(months_from_register)), months_from_register := mean(data$months_from_register, na.rm = T)]
data[which(is.na(k_var_days_between_visits_15d)), k_var_days_between_visits_15d := mean(data$k_var_days_between_visits_15d, na.rm = T)]

## clean gender categories
data[gender == "Ð–", gender := 'cat1']
data[gender == "Ðœ", gender := 'cat2']
data[gender == "", gender := 'cat3']
data[gender == "Ð\u009dÐµ Ð¾Ð¿Ñ€ÐµÐ´ÐµÐ»ÐµÐ½", gender := 'cat4']

prob_W <- 0.75

## remove missing columns
idx_col <- which(apply(data, 2,anyNA))
data <- data[, -idx_col, with = F]

## split in training and test set
idx_base <- sample(1:nrow(data), 50000, replace = F)
dt_train <- data[idx_base,]
dt_train_small <- copy(dt_train[1:round(nrow(dt_train)*0.1),])
dt_train_cal <- dt_train[1:round(nrow(dt_train)*0.8),]
dt_val <- dt_train[(round(nrow(dt_train)*0.8)+1):round(nrow(dt_train)*0.9),]
dt_cal <- dt_train[-(1:round(nrow(dt_train)*0.9)),]
dt_test <- data[-idx_base,]

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
saveRDS(results_rf_rest, paste0('results/sim_sampling_parameters/', 'lenta', '/param_rf_rest.RDS'))
saveRDS(results_rf_platt, paste0('results/sim_sampling_parameters/', 'lenta', '/param_rf_platt.RDS'))
saveRDS(results_rf_isotonic, paste0('results/sim_sampling_parameters/', 'lenta', '/param_rf_isotonic.RDS'))
saveRDS(results_rf_tuned, paste0('results/sim_sampling_parameters/', 'lenta', '/param_rf_tuned.RDS'))
saveRDS(results_rf_small, paste0('results/sim_sampling_parameters/', 'lenta', '/param_rf_small.RDS'))
fwrite(dt_test, paste0('results/sim_sampling_parameters/', 'lenta', '/dt_exp.csv'))
