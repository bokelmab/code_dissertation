if(!require("forecast")) install.packages("forecast"); library("forecast")

#wrapper function to plot two standardized time series in one graphic
plot_multi_ts = function(p_ts1,p_ts2,p_main,p_ylab=''){
  scaled1 = ts(scale(p_ts1),start=start(p_ts1),frequency=frequency(p_ts1))
  scaled2 = ts(scale(p_ts2),start=start(p_ts2),frequency=frequency(p_ts2))
  all_data <- c(scaled1,scaled2)
  plot(scaled1, main=p_main,xlab='',ylab=p_ylab,ylim=c(min(all_data),max(all_data)))
  lines(window(scaled2,start=start(scaled1)),col='red')
}


#input: time series
#output: innovations and model
prewhiten <- function(p_ts){
  model <- auto.arima(p_ts,D=1)
  return(list(innovations = model$residuals, model = model))
}

#input: two time series
#output: pre-whitened cross correlation. Filter is taken from p_ts1
prewhitened_cc <- function(p_ts1,p_ts2){
  prew1 <- prewhiten(p_ts1)
  filtered1 <- prew1[['innovations']]
  filtered2 <- Arima(p_ts2,model=prew1[['model']])$residuals
  ccf(filtered1,filtered2)
}




