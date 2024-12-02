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


get_sampling_parameters <- function(p_data, p_model, p_outcome_func = NULL){
  
  ## calculate treatment proportion and sample size
  propW <- mean(p_data$W)
  N_sample <- nrow(p_data)
  
  ## sub functions
  calc_RH <- function(p_pH, p_QV){
    return(pmin(1/(p_pH + (1-p_pH)/sqrt(p_QV)), 0.5/p_pH))
  }
  calc_var_ate_hs <- function(p_pH = pH, p_VpH, p_RH, p_VpL, p_Nsample = N_sample){
    return((p_pH*p_VpH/p_RH+((1-p_pH)^2)*p_VpL/(1-p_RH*p_pH))/p_Nsample)
  }
  calc_var_strat <- function(p_propY, p_idx_H0 = idx_H0, p_idx_H1 = idx_H1, p_idx_L0 = idx_L0,
                             p_idx_L1 = idx_L1, p_propW = propW){
    
    varH0 <- mean(p_propY[idx_H0])*mean(1-p_propY[idx_H0])
    varH1 <- mean(p_propY[idx_H1])*mean(1-p_propY[idx_H1])
    VpH <- varH0/p_propW+varH1/(1-p_propW)
    varL0 <- mean(p_propY[idx_L0])*mean(1-p_propY[idx_L0])
    varL1 <- mean(p_propY[idx_L1])*mean(1-p_propY[idx_L1])
    VpL <- varL0/p_propW+varL1/(1-p_propW)
    
    return(list(VpH = VpH, VpL = VpL))
    
  }
    
  ##### Variance of ATE estimator by random sampling ############################################
  if('mux' %in% names(p_data)){
    ## true variance 
    var_rand <- mean(p_data[W==1,]$expy)*mean(1-p_data[W==1,]$expy)/(propW*N_sample) + mean(p_data[W==0,]$expy)*mean(1-p_data[W==0,]$expy)/((1-propW)*N_sample)
    
    ## true variance, under assumption of no treatment effect 
    var_rand_ass <- mean(p_data[W==1,]$mux)*mean(1-p_data[W==1,]$mux)/(propW*N_sample)+ mean(p_data[W==0,]$mux)*mean(1-p_data[W==0,]$mux)/((1-propW)*N_sample)
  }else{
    var_rand <- NULL
    var_rand_ass <- NULL
  }
  
  ## estimated variance 
  var_rand_est <- mean(p_data[W==1,]$pred_prob)*mean(1-p_data[W==1,]$pred_prob)/(propW*N_sample) + mean(p_data[W==0,]$pred_prob)*mean(1-p_data[W==0,]$pred_prob)/((1-propW)*N_sample)
  
  ##### Calculation of Vp_H and Vp_L for different values pH
  ## for loop over pH
  pH <- c()
  Vp_L_values <- c()
  Vp_H_values <- c()
  Vp_L_values_ass <- c()
  Vp_H_values_ass <- c()
  Vp_L_values_est <- c()
  Vp_H_values_est <- c()
  for(i_p in 1:99){
    
    pH_i <- i_p/100
    pH <- c(pH, pH_i)
    
    ## separation into strata and treatment groups
    idx_H0 <- which(p_data$pred_prob > quantile(p_data$pred_prob,1-pH_i) & p_data$W == 0)
    idx_H1 <- which(p_data$pred_prob > quantile(p_data$pred_prob,1-pH_i) & p_data$W == 1)
    idx_L0 <- which(p_data$pred_prob <= quantile(p_data$pred_prob,1-pH_i) & p_data$W == 0)
    idx_L1 <- which(p_data$pred_prob <= quantile(p_data$pred_prob,1-pH_i) & p_data$W == 1)
    
    ## true variance quotient
    if('mux' %in% names(p_data)){
      erg_true <- calc_var_strat(p_propY = p_data$expy)
      Vp_L_values <- c(Vp_L_values, erg_true$VpL)
      Vp_H_values <- c(Vp_H_values, erg_true$VpH)
    }else{
      Vp_L_values <- c(Vp_L_values, NULL)
      Vp_H_values <- c(Vp_H_values, NULL)
    }
    
    
    ## true variance quotient, under assumption of no treatment effect
    if('mux' %in% names(p_data)){
      erg_ass <- calc_var_strat(p_propY = p_data$mux)
      Vp_L_values_ass <- c(Vp_L_values_ass, erg_ass$VpL)
      Vp_H_values_ass <- c(Vp_H_values_ass, erg_ass$VpH)
    }else{
      Vp_L_values_ass <- c(Vp_L_values_ass, NULL)
      Vp_H_values_ass <- c(Vp_H_values_ass, NULL)
    }
    
    ## estimated variance quotient
    erg_est <- calc_var_strat(p_propY = p_data$pred_prob)
    Vp_L_values_est <- c(Vp_L_values_est, erg_est$VpL)
    Vp_H_values_est <- c(Vp_H_values_est, erg_est$VpH)
    
  }
  
  ############## estimation of Q_H and R_H values ######################
  if('mux' %in% names(p_data)){
    Q_V <- Vp_H_values/Vp_L_values
    Q_V_ass <- Vp_H_values_ass/Vp_L_values_ass
  }else{
    Q_V <- NULL
    Q_V_ass <- NULL
  }
  Q_V_est <- Vp_H_values_est/Vp_L_values_est
  
  if('mux' %in% names(p_data)){
    R_H <- calc_RH(p_pH = pH, p_QV = Q_V) 
    R_H_ass <- calc_RH(p_pH = pH, p_QV = Q_V_ass)
  }else{
    R_H <- NULL 
    R_H_ass <- NULL
  }
  R_H_est <- calc_RH(p_pH = pH, p_QV = Q_V_est)
  
  ############## estimation of variance of ATE estimators ##############
  if('mux' %in% names(p_data)){
    var_ATE_opt <- calc_var_ate_hs(p_VpH = Vp_H_values, p_RH = R_H, p_VpL = Vp_L_values) ## true variance with optimal R_H
    var_ATE_opt_ass <- calc_var_ate_hs(p_VpH = Vp_H_values_ass, p_RH = R_H, p_VpL = Vp_L_values_ass) ## true variance with optimal R_H, under assumption of no treatment effect
    var_ATE_est_ass <- calc_var_ate_hs(p_VpH = Vp_H_values_ass, p_RH = R_H_est, p_VpL = Vp_L_values_ass) ## true variance with optimal R_H, under assumption of no treatment effect
    var_ATE_est <- calc_var_ate_hs(p_VpH = Vp_H_values_est, p_RH = R_H_est, p_VpL = Vp_L_values_est) ## true variance with estimated R_H
  }else{
    var_ATE_opt <- NULL
    var_ATE_opt_ass <- NULL
    var_ATE_est_ass <- NULL
    var_ATE_est <- NULL
  }
  est_var_ATE_est <- calc_var_ate_hs(p_VpH = Vp_H_values_est, p_RH = R_H_est, p_VpL = Vp_L_values_est) ## estimated variance with estimated R_H
  
  ### functions to calculate ATE estimator variance for given sampling parameters #####################
  calc_true_var_ate <- function(p_pH, p_RH){
    return(calc_var_ate_hs(p_pH = p_pH, p_VpH = Vp_H_values[which(pH == p_pH)], p_RH = p_RH, p_VpL = Vp_L_values[which(pH == p_pH)])) ## estimated variance with estimated R_H
  }
  calc_est_var_ate <- function(p_pH, p_RH){
    return(calc_var_ate_hs(p_pH = p_pH, p_VpH = Vp_H_values_est[which(pH == p_pH)], p_RH = p_RH, p_VpL = Vp_L_values_est[which(pH == p_pH)])) ## estimated variance with estimated R_H
  }
  calc_ass_var_ate <- function(p_pH, p_RH){
    return(calc_var_ate_hs(p_pH = p_pH, p_VpH = Vp_H_values_ass[which(pH == p_pH)], p_RH = p_RH, p_VpL = Vp_L_values_ass[which(pH == p_pH)])) ## estimated variance with estimated R_H
  }
  
  ## function for outcome predictions
  if(is.null(p_outcome_func)){ ## no calibration
    create_outcome_func <- function(){
      outcome_func <- function(p_data){
        return(predict(p_model, data = p_data)$predictions)
      }
    }
    outcome_func <- create_outcome_func()
  }else{ ## calibration
    outcome_func <- p_outcome_func
  }
  
  erg <- list(
    pH = pH,
    pH_choose = pH[which.min(est_var_ATE_est)], ### estimation of best sampling parameters
    Q_V = Q_V,
    Q_V_ass = Q_V_ass,
    Q_V_est = Q_V_est,
    R_H = R_H,
    R_H_ass = R_H_ass,
    R_H_est = R_H_est,
    R_H_choose = R_H_est[which.min(est_var_ATE_est)], ### estimation of best sampling parameters
    var_rand = var_rand,
    var_rand_est = var_rand_est,
    calc_true_var_ate = calc_true_var_ate,
    calc_ass_var_ate = calc_ass_var_ate,
    calc_est_var_ate = calc_est_var_ate,
    Vp_L_values_ass = Vp_L_values_ass,
    Vp_H_values_ass = Vp_H_values_ass,
    Vp_L_values_est = Vp_L_values_est,
    Vp_H_values_est = Vp_H_values_est,
    outcome_model = p_model,
    outcome_func = outcome_func
  )
  
  return(erg)
  
}

