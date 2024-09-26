## calculate Qini curve
calc_CATE <- function(p_Y, p_W, p_pred, p_cut){
  y_treat <- p_Y[p_pred > quantile(p_pred, p_cut) & p_W == 1]
  y_contr <- p_Y[p_pred > quantile(p_pred, p_cut) & p_W == 0]
  qini_mean <- (mean(y_treat)-mean(y_contr))*length(y_treat)
  qini_sd <- sqrt(var(y_treat)/length(y_treat)+var(y_contr)/length(y_contr))*length(y_treat)
  return(list(qini_mean = qini_mean, qini_sd = qini_sd))
}

calc_qini_curve <- function(p_preds, p_dt_test){
  
  qini_curve <- c()
  for(i_dec in 9:0){
    
    ## Qini curve and standard deviation with original target 
    res_qini <- p_dt_test %$% calc_CATE(p_Y = Y, p_W = W, p_pred = p_preds, p_cut = i_dec/10)
    qini_curve <- c(qini_curve, res_qini$qini_mean)
    
  }
  return(c(0, qini_curve))
}

