library(MASS)
library(ranger)
library(data.table)
library(magrittr)
library(latex2exp)

## required functions
source('code/helper/functions_simulation.R')

## Function to calculate R_H from outcome variance estimates
calc_RH <- function(p_pH, p_QV){
  return(pmin(1/(p_pH + (1-p_pH)/sqrt(p_QV)), 0.5/p_pH))
}


## Simulation scenario 1
samp_param_sim1 <- readRDS(paste0('results/sim_sampling_parameters/', 'sim_scen1', '/param_rf_', 'tuned', '.RDS'))
var_rand_sim1 <- samp_param_sim1$var_rand_est
var_HS_sim1 <- samp_param_sim1$calc_est_var_ate(samp_param_sim1$pH, calc_RH(samp_param_sim1$pH, samp_param_sim1$Vp_H_values_est/samp_param_sim1$Vp_L_values_est))

## Simulation scenario 2
samp_param_sim2 <- readRDS(paste0('results/sim_sampling_parameters/', 'sim_scen2', '/param_rf_', 'tuned', '.RDS'))
var_rand_sim2 <- samp_param_sim2$var_rand_est
var_HS_sim2 <- samp_param_sim2$calc_est_var_ate(samp_param_sim2$pH, calc_RH(samp_param_sim2$pH, samp_param_sim2$Vp_H_values_est/samp_param_sim2$Vp_L_values_est))

## Simulation scenario 3
samp_param_sim3 <- readRDS(paste0('results/sim_sampling_parameters/', 'sim_scen3', '/param_rf_', 'tuned', '.RDS'))
var_rand_sim3 <- samp_param_sim3$var_rand_est
var_HS_sim3 <- samp_param_sim3$calc_est_var_ate(samp_param_sim3$pH, calc_RH(samp_param_sim3$pH, samp_param_sim3$Vp_H_values_est/samp_param_sim3$Vp_L_values_est))

## Criteo
samp_param_criteo <- readRDS(paste0('results/sim_sampling_parameters/', 'criteo', '/param_rf_', 'tuned', '.RDS'))
var_rand_criteo <- samp_param_criteo$var_rand_est
var_HS_criteo <- samp_param_criteo$calc_est_var_ate(samp_param_criteo$pH, calc_RH(samp_param_criteo$pH, samp_param_criteo$Vp_H_values_est/samp_param_criteo$Vp_L_values_est))

## Lenta
samp_param_lenta <- readRDS(paste0('results/sim_sampling_parameters/', 'lenta', '/param_rf_', 'tuned', '.RDS'))
var_rand_lenta <- samp_param_lenta$var_rand_est
var_HS_lenta <- samp_param_lenta$calc_est_var_ate(samp_param_lenta$pH, calc_RH(samp_param_lenta$pH, samp_param_lenta$Vp_H_values_est/samp_param_lenta$Vp_L_values_est))



## plot results
layout(matrix(c(1,2), ncol=1, byrow=TRUE), heights = c(0.9,0.1))
par(mar = c(3, 3, 3, 3), xpd = TRUE)

plot(samp_param_sim1$pH[1:50],((1-var_HS_sim1/var_rand_sim1)*100)[1:50], type = 'l', ylab = TeX('% variance reduction'), 
     xlab = TeX('$p_H$'), lwd = 3, main = TeX('Estimated potential of HS-sampling'),
     mgp=c(1.5,0.5,0), col = 'orange', ylim = c(-3, 70))
lines(samp_param_sim2$pH[1:50],((1-var_HS_sim2/var_rand_sim2)*100)[1:50], col = 'black', lwd = 3)
lines(samp_param_sim3$pH[1:50],((1-var_HS_sim3/var_rand_sim3)*100)[1:50], col = 'red', lwd = 3)
lines(samp_param_criteo$pH[1:50],((1-var_HS_criteo/var_rand_criteo)*100)[1:50], col = 'blue', lwd = 3)
lines(samp_param_lenta$pH[1:50],((1-var_HS_lenta/var_rand_lenta)*100)[1:50], col = 'green', lwd = 3)

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('orange', 'black', 'red', 'blue', 'green')
legend(x = "bottom",inset = 0,
       legend = c('Sim 1','Sim 2', 'Sim 3', 'Criteo', 'Lenta'), 
       col=plot_colors, lty = c(1,1), lwd=5, cex=1, horiz = TRUE, bty = 'n')

