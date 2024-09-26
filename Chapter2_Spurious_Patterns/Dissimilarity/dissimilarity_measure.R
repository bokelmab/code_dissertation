#requeired files
source('WrapperFunctions/time_calculations.R')

euclidean_distance_from_mean <- function(p_series){
  
  euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  
  mean_series <- apply(p_series,1,mean)
  eucl_dist <- c()
  for(i in 1:ncol(p_series)){
    eucl_dist[i] <- euc_dist(p_series[,i],mean_series)
  }
  return(mean(eucl_dist))
}

between_ts_variance <- function(p_series,p_start,p_end){
  
  #consider specific window
  new_col <- ts(scale(window(p_series[,1],start=p_start,end=p_end)),start=p_start,frequency=12)
  ts_window <- data.frame(new_col)
  for(i in 2:ncol(p_series)){
    new_col <- ts(scale(window(p_series[,i],start=p_start,end=p_end)),start=p_start,frequency=12)
    ts_window <- data.frame(ts_window,new_col)
  }
  
  #calculate variance
  window_variance <- euclidean_distance_from_mean(ts_window)
  
  return(window_variance)
  
}

dissimilarity <- function(p_series,p_boundary,deseasonalise=F){
  
  if(deseasonalise){
    for(i in 1:ncol(p_series)){
      p_series[,i] <- p_series[,i] - stl(p_series[,i], s.window='periodic', s.degree = 1)$time.series[,'seasonal']
      p_series[,i] <- ts(scale(p_series[,i]),start=start(p_series[,1]),frequency=12)
    }
  }
  
  variance_ts <- c()
  
  for(i in (p_boundary):(nrow(p_series)-p_boundary-1)){
    variance_ts[i-p_boundary+1] <- between_ts_variance(p_series,increase_month(start(p_series[,1]),i-p_boundary),
                        increase_month(start(p_series[,1]),i+p_boundary))
  }
  
  #convert to time series
  variance_ts <- ts(variance_ts,start=increase_month(start(p_series[,1]),p_boundary),frequency=12)
  
  return(variance_ts)
}

