library(MASS)
library(ranger)
library(data.table)
library(magrittr)
library(latex2exp)

## required functions
source('code/helper/functions_simulation.R')

## plot results
layout(matrix(c(1,2), ncol=1, byrow=TRUE), heights = c(0.9,0.1))

par(mar = c(3, 3, 3, 3), xpd = TRUE)
## heteroskedasticity
pH <- 0.1
VH <- 6
VL <- 1
p <- (100:500)/1000
up_bound <- 1/(pH+(1-pH)*VL/VH)
in_green <- 1+(up_bound-1)/100*(0:100)
in_red <- up_bound + (max(p/pH)-up_bound)/100*(0:100)
default_variance <- pH*VH+(1-pH)*VL
opt_osratio <- sqrt(VH)/(pH*sqrt(VH)+(1-pH)*sqrt(VL))
est_osratio <- sqrt(VH*1.9)/(pH*sqrt(VH*1.9)+(1-pH)*sqrt(VL*0.6))



plot(p/pH,(pH^2)*VH/p+((1-pH)^2)*VL/(1-p), xaxt = "n", xlab = TeX(''), type = 'l', lwd = 3, ylab = TeX('$Var[\\hat{ATE}_{HS}]$'), 
     main = TeX('Variance of ATE estimate based on HS-sampling'), mgp=c(1,1,0), ylim = c(1,1.8), yaxt = 'n')
polygon(x = c(in_green, rev(in_green)), y = c((pH^2)*VH/in_green/pH+((1-pH)^2)*VL/(1-in_green*pH), rep(default_variance, length(in_green))), col = "#65BFFF")
polygon(x = c(in_red, rev(in_red)), y = c((pH^2)*VH/in_red/pH+((1-pH)^2)*VL/(1-in_red*pH), rep(default_variance, length(in_red))), col = "red")
lines(p/pH,(pH^2)*VH/p+((1-pH)^2)*VL/(1-p), lwd = 3)
lines(p/pH,(pH^2)*(VH*1.9)/p+((1-pH)^2)*(VL*0.6)/(1-p), lty = 1, lwd = 3, col = 'orange')
lines(c(1,1),c(0.95,default_variance), lty = 3, col = "#65BFFF")
axis(1, at=1, labels=TeX('1'), mgp=c(3,1.5,0))
lines(c(up_bound, up_bound),c(0.95,default_variance), lty = 3, col = "#65BFFF")
axis(1, at=up_bound, labels=TeX('$\\frac{1}{p_H+(1-p_H)\\cdot Q_{V}^{-1}}$'), mgp=c(3,2.5,0))
lines(c(opt_osratio,opt_osratio), c(0.95, (pH*sqrt(VH)+(1-pH)*sqrt(VL))^2), lty = 3)
axis(1, at=opt_osratio, labels=TeX('$R_H(Q_V)$'), mgp=c(3,1.5,0))
lines(c(est_osratio,est_osratio), c(0.95, (pH*sqrt(VH*1.9)+(1-pH)*sqrt(VL*0.6))^2), lty = 3)
axis(1, at=est_osratio, labels=TeX('$R_H(\\hat{Q}_V)$'), mgp=c(3,1.5,0))
arrows(opt_osratio, 1, est_osratio, 1)

##legend
par(mar=c(0, 0, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('black', 'orange')
legend(x = "bottom",inset = 0,
       legend = c(TeX('actual variance'),TeX('estimated variance')), 
       col=plot_colors, lty = c(1,1), lwd=c(3,3), cex=1, horiz = TRUE, bty = "n")

############ Graphik fuer Verteidigung #################################################
plot(p/pH,(pH^2)*VH/p+((1-pH)^2)*VL/(1-p), xlab = TeX('$R_H$'), xaxt = 'n', type = 'l', lwd = 3, ylab = TeX('$Var[\\hat{ATE}_{HS}]$'), 
     main = TeX('Variance of ATE estimate based on HS-sampling'), mgp=c(1,1,0), ylim = c(1,1.8), yaxt = 'n', cex.lab = 1.2)
polygon(x = c(in_green, rev(in_green)), y = c((pH^2)*VH/in_green/pH+((1-pH)^2)*VL/(1-in_green*pH), rep(default_variance, length(in_green))), col = "#65BFFF")
polygon(x = c(in_red, rev(in_red)), y = c((pH^2)*VH/in_red/pH+((1-pH)^2)*VL/(1-in_red*pH), rep(default_variance, length(in_red))), col = "red")
lines(p/pH,(pH^2)*VH/p+((1-pH)^2)*VL/(1-p), lwd = 5)
lines(p/pH, rep(1.5, length(p)), lwd = 5, lty = 3)
legend('bottomright', legend = c(TeX('Proportional Sampling ($R_H=1$)'), 'HS Sampling'), lty = c(3,1), lwd = c(5,5))
axis(1, at = (p/pH)[1], labels = '1')



############ Graphic variance reduction potential assessment ###########################

rel1 <- expression(X[,1]+0.5*X[,2]+X[,3]*X[,4]-4)
treat1 <- expression(0.1)
rel2 <- expression(X[,1]^2+0.5*X[,2]+X[,3]*X[,4]-7)
treat2 <- expression(0.1+(1+X[,5]))
rel3 <- expression(0.1*exp(X[,1])+0.5*X[,2]^3+X[,3]-7)
treat3 <- expression(0.1+X[,5]*X[,6])

rel <- rel3
treat <- treat3
data <- generate_data(p_n = 100000, p_rel = rel, p_treat = treat)
rf_base <- ranger(Y ~ ., data = data[W==0, -c('W'), with = F], max.depth = 5)
data$pred_prob <- predict(rf_base, data = data)$predictions
opt_split_sim3 <- find_split(data)


rel <- rel2
treat <- treat2
data <- generate_data(p_n = 100000, p_rel = rel, p_treat = treat)
rf_base <- ranger(Y ~ ., data = data[W==0, -c('W'), with = F], max.depth = 5)
data$pred_prob <- predict(rf_base, data = data)$predictions
opt_split_sim2 <- find_split(data)


rel <- rel1
treat <- treat1
data <- generate_data(p_n = 100000, p_rel = rel, p_treat = treat)
rf_base <- ranger(Y ~ ., data = data[W==0, -c('W'), with = F], max.depth = 5)
data$pred_prob <- predict(rf_base, data = data)$predictions
opt_split_sim1 <- find_split(data)


## Criteo data
data <- fread('data/criteo-uplift-v2.1.csv')
names(data)[which(names(data)== 'treatment')] <- 'W'
target <- 'conversion'
y_data <- data[, target, with = F]
data$Y <- data[, target, with = F]
data <- data[,-c('conversion', 'visit', 'exposure'), with = F]
data$Y <- y_data
idx_base <- sample(1:nrow(data), round(0.01*nrow(data)), replace = F)
dt_base <- data[idx_base,]
dt_exp <- data[-idx_base,]
rf_base <- ranger(Y ~ ., data = dt_base[W==0,-c('W'), with = F], max.depth = 5, num.trees = 500)
check <- copy(dt_exp[sample(1:nrow(dt_exp), 1000000),])
check$pred_prob <- predict(rf_base, check)$predictions
opt_split_criteo <- find_split(check)


## lenta
dt_exp <- fread('results/lenta/dt_exp.csv')
opt_split_lenta <- find_split(dt_exp)


## plot results
layout(matrix(c(1,2), ncol=1, byrow=TRUE), heights = c(0.9,0.1))
par(mar = c(3, 3, 3, 3), xpd = TRUE)

plot((0:100)/100,rev(1-opt_split_sim1$var_rel_opt)*100, type = 'l', ylab = TeX('% variance reduction'), 
     xlab = TeX('$p_H$'), lwd = 3, main = TeX('Estimated potential of HS-sampling'),
     mgp=c(1.5,0.5,0), col = 'orange', ylim = c(0, 66))
lines((0:100)/100,rev(1-opt_split_sim2$var_rel_opt)*100, lwd = 3)
lines((0:100)/100,rev(1-opt_split_sim3$var_rel_opt)*100, col = 'red', lwd = 3)
lines((0:100)/100,rev(1-opt_split_criteo$var_rel_opt)*100, col = 'blue', lwd = 3)
lines((0:100)/100,rev(1-opt_split_lenta$var_rel_opt)*100, col = 'green', lwd = 3)

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('orange', 'black', 'red', 'blue', 'green')
legend(x = "bottom",inset = 0,
       legend = c('Sim 1','Sim 2', 'Sim 3', 'Criteo', 'Lenta'), 
       col=plot_colors, lty = c(1,1), lwd=5, cex=1, horiz = TRUE, bty = 'n')

