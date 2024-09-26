gt_analysis <- function(p_gt,p_arrival){
  
  #pre-whitened cross-correlation
  results_pcc <- prewhitened_cc(p_gt,p_arrival)
  cc_lags <- results_pcc$acf[which(results_pcc$lag>=-0.5 & results_pcc$lag<=0)] #cc till lag 6
  
  #inclusion decision 
  cut_off <- 2 / sqrt(length(p_arrival)) #2 times standard deviation of sample ccf if no cross-correlation
  if(cc_lags[7] > cut_off){
    include <- 1
  }else{
    include <- 0
  }
  
  #lag order selection
  lags <- 7 - which(cc_lags > cut_off)
  
  return(list(include=include,lags=lags,cross_correlations=results_pcc$acf))
}