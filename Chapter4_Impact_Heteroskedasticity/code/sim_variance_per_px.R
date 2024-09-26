library(data.table)
library(dplyr)
library(magrittr)
library(grf)
library(ranger)
library(latex2exp)


## Ab hier relevant fuer Paper
dist_preds_cf <- NULL
dist_preds_t <- NULL
dist_preds_s <- NULL
for(i_sim in 1:10000){
  
  x_train <- (0:1000)/1000
  y0_train <- rbinom(n = length(x_train), size = 1, prob = x_train)
  y1_train <- rbinom(n = length(x_train), size = 1, prob = x_train)
  dt_train <- data.table(y = c(y0_train, y1_train), x = c(x_train, x_train), w = c(rep(0, 1001), rep(1, 1001)))
  cf_model <- causal_forest(X = dt_train[,-c('w', 'y'), with = F], Y = dt_train$y, W = dt_train$w, W.hat = 0.5, num.trees = 500)
  rf0 <- ranger(y ~ x, data = dt_train[w == 0,])
  rf1 <- ranger(y ~ x, data = dt_train[w == 1,])
  rfs <- ranger(y ~ x + w, data = dt_train)
  preds_cf <- predict(cf_model, data.table(x = (0:100)/100))$predictions
  preds_t <- predict(rf1, data = data.table(x = (0:100)/100))$predictions - predict(rf0, data = data.table(x = (0:100)/100))$predictions
  preds_s <- predict(rfs, data = data.table(x = (0:100)/100, w = 1))$predictions - predict(rfs, data = data.table(x = (0:100)/100, w = 0))$predictions
  
  dist_preds_cf %<>% rbind(preds_cf)
  dist_preds_t %<>% rbind(preds_t)
  dist_preds_s %<>% rbind(preds_s)
  if(i_sim %% 1000 == 0){
    print(i_sim)
  }
}

## save results
saveRDS(dist_preds_cf, 'results/dist_preds_cf.RDS')
saveRDS(dist_preds_t, 'results/dist_preds_t.RDS')
saveRDS(dist_preds_s, 'results/dist_preds_s.RDS')

dist_preds_cf <- readRDS('results/dist_preds_cf.RDS')
dist_preds_s <- readRDS('results/dist_preds_s.RDS')
dist_preds_t <- readRDS('results/dist_preds_t.RDS')

#################### Graphic paper ###############################

######### Calculate results ####################################################

## Ranking
order_tepreds_cf <- apply(-dist_preds_cf, 1, order)
order_tepreds_t <- apply(-dist_preds_t, 1, order)
order_tepreds_s <- apply(-dist_preds_s, 1, order)

## specify graphic settings
layout(matrix(c(1,2,3,4), ncol=1, byrow=T), heights = c(0.1, 0.23, 0.23, 0.23))
par(mar=c(3, 2.5, 0.5, 2), mgp=c(1.3, 0.4, 0))

plot.new()
text(0.5,0.5,"Distribution of CATE estimates",cex=2,font=2)

## CRF
exp_values <- apply(dist_preds_cf, 2, mean)
var_values <- apply(dist_preds_cf, 2, var)
plot((0:100)/100, exp_values, ylim = c(min(exp_values), max(sqrt(var_values))), type = 'l', 
     lwd = 3, xlab = TeX('$x$'), ylab = '',col = 'red')
lines((0:100)/100, sqrt(var_values), lty = 1, lwd = 3, col = 'blue')
text(median((0:100)/100), 0.05 , 'CRF', cex = 1.5)

## s-learner
exp_values <- apply(dist_preds_s, 2, mean)
var_values <- apply(dist_preds_s, 2, var)
plot((0:100)/100, exp_values, ylim = c(min(exp_values), max(sqrt(var_values))), type = 'l', 
     lwd = 3, xlab = TeX('$x$'), ylab = '',col = 'red')
lines((0:100)/100, sqrt(var_values), lty = 1, lwd = 3, col = 'blue')
text(median((0:100)/100), 0.03 , 'S-learner', cex = 1.5)

## t-learner
exp_values <- apply(dist_preds_t, 2, mean)
var_values <- apply(dist_preds_t, 2, var)
plot((0:100)/100, exp_values, ylim = c(min(exp_values), max(sqrt(var_values))), type = 'l', 
     lwd = 3, xlab = TeX('$x$'), ylab = '',col = 'red')
lines((0:100)/100, sqrt(var_values), lty = 1, lwd = 3, col = 'blue')
text(median((0:100)/100), 0.2 , 'T-learner', cex = 1.5)

## plot u-shaped rf percentiles per percentile of ranking

#scatter.smooth((0:50)/100, (apply(order_tepreds_x, 1, mean)-1)/100, span = 1/3, lwd = 1, lpars = list(col = "blue", lwd = 3, lty = 1), xaxt='n', yaxt='n')
layout(matrix(c(1,2,3,4), ncol=1, byrow=T), heights = c(0.1, 0.23, 0.23, 0.23))
par(mar=c(3, 2.5, 0.5, 2), mgp=c(1.3, 0.4, 0))

plot.new()
text(0.5,0.5,"Ranking",cex=2,font=2)
av_var_quant_cf <- apply(((order_tepreds_cf-1)/100)*(1-(order_tepreds_cf-1)/100),1,mean)  ## expected varaince at ranking position
plot((0:100)/100, av_var_quant_cf, col = "black", lwd = 3, lty = 1, xlab = TeX('$F_{\\hat{\\tau}}$'), ylab = TeX('$E[\\mu_{x}(1-\\mu_{x})|F_{\\hat{\\tau}}]$'), type = 'l')
text(median((0:100)/100), 0.17, 'CRF', cex = 1.5)
av_var_quant_s <- apply(((order_tepreds_s-1)/100)*(1-(order_tepreds_s-1)/100),1,mean)  ## expected varaince at ranking position
plot((0:100)/100, av_var_quant_s, col = "black", lwd = 3, lty = 1, xlab = TeX('$F_{\\hat{\\tau}}$'), ylab = TeX('$E[\\mu_{x}(1-\\mu_{x})|F_{\\hat{\\tau}}]$'), type = 'l')
text(median((0:100)/100), 0.17 , 'S-learner', cex = 1.5)
av_var_quant_t <- apply(((order_tepreds_t-1)/100)*(1-(order_tepreds_t-1)/100),1,mean)  ## expected varaince at ranking position
plot((0:100)/100, av_var_quant_t, col = "black", lwd = 3, lty = 1, xlab = TeX('$F_{\\hat{\\tau}}$'), ylab = TeX('$E[\\mu_{x}(1-\\mu_{x})|F_{\\hat{\\tau}}]$'), type = 'l')
text(median((0:100)/100), 0.14 , 'T-learner', cex = 1.5)

plot.new()


