##required libraries
library(data.table)
library(dplyr)
library(magrittr)
library(ranger)
library(grf)
library(latex2exp)

## required helper
source('code/helper/Qini_calculation.R')

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
data <- as.data.table(model.matrix(Y~-1+., data))
data$Y <- y_data

## split in training and test set
idx_train <- sample(1:nrow(data), round(0.1*nrow(data)), replace = F)
dt_train <- data[idx_train,]
dt_test <- data[-idx_train,]

## build random forest models for control and intervention group, for all values of min.node.size
rf0_10 <- ranger(Y ~ ., data = dt_train[W==0,-c('W'), with = F], min.node.size = 10, num.trees = 500)
rf0_100 <- ranger(Y ~ ., data = dt_train[W==0,-c('W'), with = F], min.node.size = 100, num.trees = 500)
rf0_1000 <- ranger(Y ~ ., data = dt_train[W==0,-c('W'), with = F], min.node.size = 1000, num.trees = 500)

rf1_10 <- ranger(Y ~ ., data = dt_train[W==1,-c('W'), with = F], min.node.size = 10, num.trees = 500)
rf1_100 <- ranger(Y ~ ., data = dt_train[W==1,-c('W'), with = F], min.node.size = 100, num.trees = 500)
rf1_1000 <- ranger(Y ~ ., data = dt_train[W==1,-c('W'), with = F], min.node.size = 1000, num.trees = 500)

## make predictions on the test data (split data in 10 parts for computational reasons)
pred_rf0_10 <- c()
pred_rf0_100 <- c()
pred_rf0_1000 <- c()
pred_rf1_10 <- c()
pred_rf1_100 <- c()
pred_rf1_1000 <- c()
for(i_fold in 1:12){
  pred_rf0_10 <- c(pred_rf0_10, predict(rf0_10, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions)
  pred_rf0_100 <- c(pred_rf0_100, predict(rf0_100, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions)
  pred_rf0_1000 <- c(pred_rf0_1000, predict(rf0_1000, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions)
  pred_rf1_10 <- c(pred_rf1_10, predict(rf1_10, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions)
  pred_rf1_100 <- c(pred_rf1_100, predict(rf1_100, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions)
  pred_rf1_1000 <- c(pred_rf1_1000, predict(rf1_1000, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions)
  
}
pred_rf0_10 <- c(pred_rf0_10, predict(rf0_10, dt_test[12000001:nrow(dt_test),])$predictions)
pred_rf0_100 <- c(pred_rf0_100, predict(rf0_100, dt_test[12000001:nrow(dt_test),])$predictions)
pred_rf0_1000 <- c(pred_rf0_1000, predict(rf0_1000, dt_test[12000001:nrow(dt_test),])$predictions)
pred_rf1_10 <- c(pred_rf1_10, predict(rf1_10, dt_test[12000001:nrow(dt_test),])$predictions)
pred_rf1_100 <- c(pred_rf1_100, predict(rf1_100, dt_test[12000001:nrow(dt_test),])$predictions)
pred_rf1_1000 <- c(pred_rf1_1000, predict(rf1_1000, dt_test[12000001:nrow(dt_test),])$predictions)

## save predictions
file_path <- paste0('simulation_results/criteo')
fwrite(data.table(pred_rf0_10, pred_rf0_100, pred_rf0_1000, pred_rf1_10, pred_rf1_100, pred_rf1_1000), paste0(file_path, '/preds_test.csv'))
fwrite(dt_test, paste0(file_path, '/dt_test.csv'))



#rm(pred_rf, pred_rf0, pred_rf1, pred_cf, pred_rfs, pred_x, pred_t, rf, rf_s, rf0, rf1, rfX0, rfx1, cf)


preds_test <- fread('simulation_results/criteo/node_size1000/preds_test.csv')
dt_test <- fread('simulation_results/criteo/node_size1000/dt_test.csv')

######################## Tuning #######################################

### best tuning parameters were node_size=1000
## NOTE: Die Laenge des Datensatzes fuer die Praediktionen fuer verschiedene node sizes
##       stimmte in den abgespeciherten Daten nicht mit dem Test-Set überein.

# idx_val <- sample(1:nrow(dt_test), length(idx_train), replace = F)
# dt_val <- dt_test[idx_val,]
# dt_test <- dt_test[-idx_val,]
# preds_val <- preds_test[idx_val,]
# preds_test <- preds_test[-idx_val,]
# 
# calc_mse <- function(p_outcome, p_preds){
#   return(mean((p_outcome-p_preds)^2))
# }
# mse_res0 <- c('size_10' = calc_mse(dt_val[W==0,]$Y, preds_val[dt_val$W == 0,]$pred_rf0_10), 
#               'size_100' = calc_mse(dt_val[W==0,]$Y, preds_val[dt_val$W == 0,]$pred_rf0_100), 
#               'size_1000' = calc_mse(dt_val[W==0,]$Y, preds_val[dt_val$W == 0,]$pred_rf0_1000)) 
# mse_res1 <- c('size_10' = calc_mse(dt_val[W==1,]$Y, preds_val[dt_val$W == 1,]$pred_rf1_10), 
#               'size_100' = calc_mse(dt_val[W==1,]$Y, preds_val[dt_val$W == 1,]$pred_rf1_100), 
#               'size_1000' = calc_mse(dt_val[W==1,]$Y, preds_val[dt_val$W == 1,]$pred_rf1_1000))



######################### Graphics paper ##############################
dt_test$pred_rf0 <- preds_test$pred_rf0
dt_test$pred_rf1 <- preds_test$pred_rf1
dt_test$pred_t <- preds_test$pred_rf1 - preds_test$pred_rf0

## get percentile of rf prediction
dt_test$perc_rf0 <- 1
for(i in 1:99){
  dt_test[pred_rf0 > quantile(dt_test$pred_rf0, i*0.01),]$perc_rf0 <- (i+1)
}

## get percentile of T-learner prediction
dt_test$perc_t <- 1
for(i in 1:99){
  dt_test[pred_t > quantile(dt_test$pred_t, i*0.01),]$perc_t <- (i+1)
}


## claculate Qini curve
calc_cum_effect <- function(p_perc_sel, p_perc_model){
  data_perc <- dt_test[which(p_perc_model >= p_perc_sel),]
  return((mean(data_perc[W==1,]$Y)-mean(data_perc[W==0,]$Y))*(nrow(data_perc)))
}
qini_t <- sapply(100:1, calc_cum_effect, p_perc_model = dt_test$perc_t)
qini_rf0 <- sapply(100:1, calc_cum_effect, p_perc_model = dt_test$perc_rf0)

plot((0:100)/100, c(0, qini_t), type = 'l', lwd = 3, col = 'red', ylab = 'cumulative incremental gains', main = 'Qini curve', xlab = TeX('$1-F_{\\hat{\\tau}}$'))
lines((0:100)/100, c(0, qini_rf0), col = 'blue', lwd = 3)
lines(c(0,1), c(0,max(qini_t)), lty = 3, lwd = 3)
legend('bottomright', legend = c('T-learner', 'Outcome model'), col = c('red', 'blue'), lty = 1)

## calculate 5% and 95% quantiles within each rf percentile
summary_t <- dt_test %>% group_by(perc_rf0) %>% summarise(mean = mean(pred_t), q95 = quantile(pred_t, 0.95), q05 = quantile(pred_t, 0.05))

## specify graphic settings
layout(matrix(c(1,2,3,4,4,4), ncol=3, byrow=T), heights = c(0.85,0.1))
par(mar=c(3, 5, 3, 3), mgp=c(1.5, 0.4, 0))


## treatment effect per response rate
#resp0 <- dt_test[W==0,] %>% group_by(perc_rf0) %>% summarise(mean_resp = mean(Y))
#resp1 <- dt_test[W==1,] %>% group_by(perc_rf0) %>% summarise(mean_resp = mean(Y))
plot(resp1$perc_rf0/100,resp1$mean_resp-resp0$mean_resp, main = 'Average treatment effect', type = 'l', lwd = 3,
     xlab = TeX('$F_{\\hat{\\mu}}$'), ylab = TeX('$\\hat{E}[\\tau_{x}|F_{\\hat{\\mu}}]$'), col = 'black', cex.main = 2, cex.lab = 2)

## noise per response rate
#pred0 <- dt_test[W==0,] %>% group_by(perc_rf0) %>% summarise(mean_resp = mean(pred_rf0))
plot(pred0$perc_rf0/100,pred0$mean_resp*(1-pred0$mean_resp), main = 'Unexplained variance', type = 'l', lwd = 3,
     xlab = TeX('$F_{\\hat{\\mu}}$'), ylab = TeX('$\\hat{E}[\\mu_{x}(1-\\mu_{x})|F_{\\hat{\\mu}}]$'), col = 'black', cex.main = 2, cex.lab = 2)

## distribution of T-learner predictions
plot(summary_t$perc_rf0/100, summary_t$mean, col = 'black', type = 'l', lwd = 3, 
     ylim = c(min(summary_t$q05), max(summary_t$q95)), xlab = TeX('$F_{\\hat{\\mu}}$'), ylab = TeX('\\hat{\\tau}'), 
     main = 'Distribution of CATE estimates', cex.main = 2, cex.lab = 2)
lines(summary_t$perc_rf0/100, summary_t$q95, lty = 3, lwd = 3, col = 'green')
lines(summary_t$perc_rf0/100, summary_t$q05, lty = 3, lwd = 3, col = 'green')

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('black', 'green')
legend(x = "bottom",inset = 0,
       legend = c('Conditional expected value', '5%- and 95%-quantiles'), 
       lty = c(1,3), lwd=3, cex=2, horiz = TRUE, seg.len = 5, text.width = 0.2, bty = 'n',
       col = plot_colors)
plot.new()

