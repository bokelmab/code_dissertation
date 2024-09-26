#input: date of format c(year, month)
#output: original date increased by p_nr_month month
increase_month = function(p_date,p_nr_month=1){
  month = p_date[2]
  year = p_date[1]
  year = year + (month-1+p_nr_month)%/%12
  month = (month+p_nr_month)%%12
  if(month==0){
    month = 12
  }
  return(c(year,month))
}

#input: two dates
#output: length of interval between dates in months
time_between <- function(p_start,p_end){
  return((p_end-p_start)[1] * 12 + (p_end-p_start)[2] + 1)
}

#input: a value, start point and end point
#output: constant time series between start and end point
constant_ts <- function(p_start,p_end,p_value){
  
  return(ts(rep(p_value,time_between(p_start,p_end)),start=p_start,end=p_end,frequency=12))
}

#input: two adjacent time series
#output: input series merged to one series
merge_ts <- function(p_ts1,p_ts2){
  if(!(increase_month(end(p_ts1),1)[1]==start(p_ts2)[1] & increase_month(end(p_ts1),1)[2]==start(p_ts2)[2])){
    print('Error: Time series are not adjacent.')
    return()
  }else{
    return(ts(c(p_ts1,p_ts2),start=start(p_ts1),frequency=12))
  }
}

#input: time of format c(year,month)
#output: time of format time(time series)
convert_time <- function(p_time){
  return(p_time[1]+(p_time[2]-1)/12)
}

#input: time series
#output: scaled time series
scale_ts <- function(p_ts){
  return(ts(scale(p_ts),start=start(p_ts),frequency=12))
}

#wrapper function to find the right index in analyse interval. Necessary because there
#seems to be a bug in R time series handling
ts_index <- function(p_time_vector,p_time){
  return(max(which(p_time_vector<convert_time(p_time))))
}

mark_point <- function(p_point,p_col='black',p_lty=2){
  abline(v=convert_time(p_point),lty=p_lty,col=p_col)
}

#construct data frame with dependend variable and lagged predictors
lagged_predictors <- function(p_arrival,p_gt,p_lags=NULL){
  
    gt <- window(p_gt,start=start(p_arrival))
    data <- data.frame(arrival=p_arrival,gt=gt)
  
  #add lagged predictors
  if(length(p_lags)>1){
    for(i in 1:(length(p_lags)-1)){
      gt_lag <- p_gt
      
      gt_lag <- stats::lag(p_gt,-p_lags[i])
      gt_lag <- window(gt_lag,start=start(p_arrival),end=end(gt))
      
      #add to data
      data <- data.frame(data,gt_lag)
      names(data)[2+i] <- paste('gt_l',p_lags[i],sep ="")
      
    }
  }
  
  
  return(data)
}

#get p_ts1 at the same time p_ts2 is observed
window_as <- function(p_ts1,p_ts2){
  return(window(p_ts1,start=start(p_ts2),end=end(p_ts2)))
}

