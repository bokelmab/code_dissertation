library(MASS)
library(ranger)
library(data.table)
library(magrittr)
library(latex2exp)
library(grf)
library(bartCause)

## required functions
source('code/helper/functions_simulation.R')

rel1 <- expression(2*(X[,1]>1)+0.5*X[,2]+X[,3]*X[,4]-4)
treat1 <- expression(X[,2]^2)
rel2 <- expression(X[,1]^2+0.5*X[,2]+X[,3]*X[,4]-7)
treat2 <- expression(0.1+(1+X[,5]))
rel3 <- expression(0.1*exp(X[,1])+0.5*X[,2]^3+X[,3]-7)
treat3 <- expression(0.1+X[,5]*X[,6])

rel <- rel3
treat <- treat3
scen <- 'sim_scen3'

get_pH_Rh <- function(p_data){
  
  ##theoretical optimum
  pH<- c()
  Rh <- c()
  var_rel_opt <- c()
  
  var_rand <- mean(p_data[W==1,]$pred_prob)*mean(1-p_data[W==1,]$pred_prob)+mean(p_data[W==0,]$pred_prob)*mean(1-p_data[W==0,]$pred_prob)
  
  ## split in high and low variance group
  for(i_p in 0:99){
    
    p <- i_p/100
    
    varH0 <- mean(p_data[pred_prob > quantile(pred_prob,p) & W==0,]$pred_prob)*mean(1-p_data[pred_prob > quantile(pred_prob,p) & W==0,]$pred_prob)
    varH1 <- mean(p_data[pred_prob > quantile(pred_prob,p) & W==1,]$pred_prob)*mean(1-p_data[pred_prob > quantile(pred_prob,p) & W==1,]$pred_prob)
    varH <- varH0+varH1
    varL0 <- mean(p_data[pred_prob <= quantile(pred_prob,p) & W==0,]$pred_prob)*mean(1-p_data[pred_prob <= quantile(pred_prob,p) & W==0,]$pred_prob)
    varL1 <- mean(p_data[pred_prob <= quantile(pred_prob,p) & W==1,]$pred_prob)*mean(1-p_data[pred_prob <= quantile(pred_prob,p) & W==1,]$pred_prob)
    varL <- varL0+varL1
    if(p == 0){
      varL <- 0
    }
    if(p==1){
      varH <- 0
    }
    
    ## optimal allocation
    pH <- c(pH, 1-p)
    
    ## proportional allocation
    Rh <- c(Rh, (1-p+p/sqrt(varH/varL))^(-1))
    
    ## variance reduction
    var_rel_opt <- c(var_rel_opt, (((1-p)*sqrt(varH)+p*sqrt(varL))^2)/var_rand)
    
  }
  return(list(pH=pH, Rh=Rh, var_red = var_rel_opt))
}


###################### Obtain groups and allocation ratio ######################


## plot results
layout(matrix(c(1,2,3,4), ncol=1, byrow=TRUE), heights = c(0.1, 0.3, 0.3, 0.3))


## title
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5, TeX("$R_H(\\hat{Q}_V)$ for a given $p_H$"),cex=2,font=1.5)

par(mar=c(3, 3, 2, 2))
data1 <- fread('results/sim_scen1/data_complete.csv')
pH_Rh1 <- get_pH_Rh(data1)
plot(pH_Rh1$pH, pH_Rh1$Rh, type = 'l', xlab = TeX('$p_H$'), ylab = TeX('$R_H(\\hat{Q}_V)$'), 
     col = 'orange', lwd = 3, main = 'Simulation 1', mgp=c(1.5,0.5,0))

par(mar=c(3, 3, 2, 2))
data2 <- fread('results/sim_scen2/data_complete.csv')
pH_Rh2 <- get_pH_Rh(data2)
plot(pH_Rh2$pH, pH_Rh2$Rh, type = 'l', xlab = TeX('$p_H$'), ylab = TeX('$R_H(\\hat{Q}_V)$'), 
     col = 'black', lwd = 3, main = 'Simulation 2', mgp=c(1.5,0.5,0))

par(mar=c(3, 3, 2, 2))
data3 <- fread('results/sim_scen3/data_complete.csv')
pH_Rh3 <- get_pH_Rh(data3)
plot(pH_Rh3$pH, pH_Rh3$Rh, type = 'l', xlab = TeX('$p_H$'), ylab = TeX('$R_H(\\hat{Q}_V)$'), 
     col = 'red', lwd = 3, main = 'Simulation 3', mgp=c(1.5,0.5,0))

