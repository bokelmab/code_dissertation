generate_data <- function(p_n, p_rel, p_treat){
  
  X <- mvrnorm(n = p_n, mu = rep(0, 10), Sigma = diag(10))
  W <- rbinom(n = p_n, size = 1, prob = 0.5)
  
  sigmoid <- function(t){
    exp(t)/(1+exp(t))
  }
  mux <- sigmoid(eval(p_rel))
  taux <- sigmoid(eval(p_rel)+eval(p_treat))-mux
  expy <- mux + W*taux
  Y <- rbinom(n = p_n, size = 1, prob = expy)
  data <- data.table(Y,W,X,mux,taux,expy)
  
  return(data)
}

find_optimal_split <- function(p_data){
  
  ##theoretical optimum
  var_rel_opt <- c()
  var_rel_prop <- c()
  ## random sampling
  var_rand <- mean(p_data[W==1,]$expy)*mean(1-p_data[W==1,]$expy)+mean(p_data[W==0,]$expy)*mean(1-p_data[W==0,]$expy)
  ## split in high and low variance group
  for(i_p in 0:100){
    
    p <- i_p/100
    
    varH0 <- mean(p_data[mux > quantile(mux,p) & W==0,]$expy)*mean(1-p_data[mux > quantile(mux,p) & W==0,]$expy)
    varH1 <- mean(p_data[mux > quantile(mux,p) & W==1,]$expy)*mean(1-p_data[mux > quantile(mux,p) & W==1,]$expy)
    varH <- varH0+varH1
    varL0 <- mean(p_data[mux <= quantile(mux,p) & W==0,]$expy)*mean(1-p_data[mux <= quantile(mux,p) & W==0,]$expy)
    varL1 <- mean(p_data[mux <= quantile(mux,p) & W==1,]$expy)*mean(1-p_data[mux <= quantile(mux,p) & W==1,]$expy)
    varL <- varL0+varL1
    if(p == 0){
      varL <- 0
    }
    if(p==1){
      varH <- 0
    }
    
    ## optimal allocation
    var_rel_opt <- c(var_rel_opt, (((1-p)*sqrt(varH)+p*sqrt(varL))^2)/var_rand)
    
    ## proportional allocation
    var_rel_prop <- c(var_rel_prop, ((1-p)*varH+p*varL)/var_rand)
    
  }
  return(list(var_rel_opt=var_rel_opt, var_rel_prop=var_rel_prop))
}

find_split <- function(p_data){
  
  ##theoretical optimum
  var_rel_opt <- c()
  var_rel_prop <- c()
  ## random sampling
  var_rand <- mean(p_data[W==1,]$pred_prob)*mean(1-p_data[W==1,]$pred_prob)+mean(p_data[W==0,]$pred_prob)*mean(1-p_data[W==0,]$pred_prob)
  ## split in high and low variance group
  for(i_p in 0:100){
    
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
    var_rel_opt <- c(var_rel_opt, (((1-p)*sqrt(varH)+p*sqrt(varL))^2)/var_rand)
    
    ## proportional allocation
    var_rel_prop <- c(var_rel_prop, ((1-p)*varH+p*varL)/var_rand)
    
  }
  return(list(var_rel_opt=var_rel_opt, var_rel_prop=var_rel_prop))
}