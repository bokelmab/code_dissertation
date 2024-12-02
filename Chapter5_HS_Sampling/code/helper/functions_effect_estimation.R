calc_ate_het <- function(p_dt_complete, p_dt_het){
  
  pH <- mean(p_dt_complete$var_group == 2)
  
  ate1 <- mean(p_dt_het[W==1 & var_group == 1,]$Y)-mean(p_dt_het[W==0 & var_group == 1,]$Y)
  ate2 <- mean(p_dt_het[W==1 & var_group == 2,]$Y)-mean(p_dt_het[W==0 & var_group == 2,]$Y)
  est_ate_het <- (1-pH)*ate1+pH*ate2
  
  var_ate1 <- var(p_dt_het[W==1 & var_group == 1,]$Y)/nrow(p_dt_het[W==1 & var_group == 1,])+var(p_dt_het[W==0 & var_group == 1,]$Y)/nrow(p_dt_het[W==0 & var_group == 1,])
  var_ate2 <- var(p_dt_het[W==1 & var_group == 2,]$Y)/nrow(p_dt_het[W==1 & var_group == 2,])+var(p_dt_het[W==0 & var_group == 2,]$Y)/nrow(p_dt_het[W==0 & var_group == 2,])
  var_ate_het <- ((1-pH)^2)*var_ate1+(pH^2)*var_ate2
  
  return(list(est_ate = est_ate_het, est_var = var_ate_het))
}

calc_auq <- function(p_preds, p_dt_test){
  ecdf_preds <- ecdf(p_preds)
  p_dt_test$dec <- ceiling(ecdf_preds(p_preds)*10)
  qini <- c()
  for(i_dec in 10:1){
    qini <- c(qini, (mean(p_dt_test[dec >= i_dec & W== 1,]$Y)-mean(p_dt_test[dec >= i_dec & W== 0,]$Y))*(11-i_dec)/10)
  }
  
  return(mean(qini))
}