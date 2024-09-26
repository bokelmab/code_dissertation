#load gt data
gt_data_germany_final <- readRDS("Data/GT_Data/germany/gt_data_germany_final.RDS")
gt_data_international_final <- readRDS("Data/GT_Data/international/gt_data_international_final.RDS")

#take average 
germany_av <- ts(apply(gt_data_germany_final,1,mean),start=c(2008,1),frequency=12)
international_av <- ts(apply(gt_data_international_final,1,mean),start=c(2008,1),frequency=12)

#generate different versions of filtered averages
germany_av_filtered2 <- stl(germany_av, s.window=7, s.degree = 1)$time.series[,'trend']
international_av_filtered2 <- stl(international_av, s.window=7, s.degree = 1)$time.series[,'trend']

#modify gt data
modify_gt_series <- function(p_series,p_divisor){
  output <- p_series
  for(i in 1:ncol(p_series)){
    output[,i] <- p_series[,i] / p_divisor * mean(p_divisor)
  }
  return(output)
}

germany_modified1 <- modify_gt_series(gt_data_germany_final,germany_av)
germany_modified2 <- modify_gt_series(gt_data_germany_final,germany_av_filtered2)

international_modified1 <- modify_gt_series(gt_data_international_final,international_av)
international_modified2 <- modify_gt_series(gt_data_international_final,international_av_filtered2)

#save data
saveRDS(germany_modified1,'Data/GT_Data/modified/germany_modified1.RDS')
saveRDS(germany_modified2,'Data/GT_Data/modified/germany_modified2.RDS')

saveRDS(international_modified1,'Data/GT_Data/modified/international_modified1.RDS')
saveRDS(international_modified2,'Data/GT_Data/modified/international_modified2.RDS')

