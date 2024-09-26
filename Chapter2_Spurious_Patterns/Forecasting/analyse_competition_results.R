#required files 
source('WrapperFunctions/cross_correlations.R')

env_competition_results <- environment()

##########Load results of for loop in file "workflow_competition"

gt_regions <- c('germany','international') 
versions <- c('','_modified1','_modified2')
list_gt_types <- c()
for(i in 1:length(gt_regions)){
  for(j in 1:length(versions)){
    list_gt_types <- c(list_gt_types,paste(gt_regions[i],versions[j],sep=''))
  }
}

#read competition results
in_sample_mape <- readRDS('Forecasting/Results/Errors/in_sample_mape.RDS')
in_sample_mase <- readRDS('Forecasting/Results/Errors/in_sample_mase.RDS')
in_sample_rmse <- readRDS('Forecasting/Results/Errors/in_sample_rmse.RDS')
out_sample_mape <- readRDS('Forecasting/Results/Errors/out_sample_mape.RDS')
out_sample_mase <- readRDS('Forecasting/Results/Errors/out_sample_mase.RDS')
out_sample_rmse <- readRDS('Forecasting/Results/Errors/out_sample_rmse.RDS')
inclusion_decision <- readRDS('Forecasting/Results/Errors/inclusion_decision.RDS')

for(i in 1:length(list_gt_types)){
  #cross correlation
  name <- paste('cross_correlations_',list_gt_types[i],sep='')
  data <- readRDS(paste('Forecasting/Results/Cross_Correlations/',name,'.RDS',sep=''))
  assign(name,data,env_competition_results)
  
  #coefficient information
  name <- paste('coefficient_information_',list_gt_types[i],sep='')
  data <- readRDS(paste('Forecasting/Results/Coefficient_Information/',name,'.RDS',sep=''))
  assign(name,data,env_competition_results) 
}


dm_results <- readRDS('Forecasting/Results/DM_Results/dm_results.RDS')

#get gt data
gt_data_germany_final <- readRDS("Data/GT_Data/germany/gt_data_germany_final.RDS")
arrival_data_matched <- readRDS("Data/Arrival_Data/arrival_data_matched.RDS")
gt_regions <- c('germany','international') 
versions <- c('','_modified1','_modified2')
gt_data_path <- 'Data/GT_Data/'
for(i in 1:length(gt_regions)){
  for(j in 2:length(versions)){
    path <- paste(gt_data_path,'modified','/',gt_regions[i],versions[j],'.RDS',sep='')
    name <- paste('gt_',gt_regions[i],'_data',versions[j],sep='')
    assign(name,readRDS(path),env_competition_results)
  }
}

##########helper functions

#comparison raw against modified 
compare_statistics <- function(p_errors,p_model1,p_model2){
  
  #means
  print('Comparison of means:')
  print(paste(p_model1,':',mean(p_errors[,p_model1]),';',p_model2,':',mean(p_errors[,p_model2])))
  
  print('Comparison of medians:')
  print(paste(p_model1,':',median(p_errors[,p_model1]),';',p_model2,':',median(p_errors[,p_model2])))
  
  #dm test
  #coding for dm test data frame
  dm_coding <- function(p_model1,p_model2){
    return(paste(p_model1,'_vs_',p_model2,sep=''))
  }
  first_better <- length(which(dm_results[,dm_coding(p_model1,p_model2)]<0.05))
  second_better <- length(which(dm_results[,dm_coding(p_model2,p_model1)]<0.05))
  print(paste(nrow(p_errors),'series;',first_better,'times was',p_model1,'superior;',second_better,'times was',p_model2,'superior'))
}

#aggregate according to region
matching_frame <- readRDS("Data/Arrival_Data/matching_frame.RDS")
match_names <- function(p_query){
  data <- get('matching_frame',env_competition_results)
  return(as.character(data[,p_query]))
}

select_queries <- function(p_in_sample,p_model){
  
  regions <- unlist(lapply(rownames(p_in_sample),match_names))
  regions_unique <- unique(regions)
  
  #get queries with minimum in-sample criterion in each region
  selected_queries <- c()
  for(i in 1:length(regions_unique)){
    in_sample_region <- p_in_sample[which(regions==regions_unique[i]),]
    selected_queries[i] <- rownames(in_sample_region)[which(in_sample_region[,p_model]==min(in_sample_region[,p_model]))]
  }
  return(selected_queries)
  
}

#function to compare models forecasting results aggregated according to region
results_region <- function(p_in_sample,p_out_sample,p_dm,p_gt_model1,p_gt_model2,p_gt_model3,p_print_dm){
  
  env_results_region <- environment()
  
  ########error results
  
  #select one query for each region
  selected_queries1 <- select_queries(p_in_sample,p_gt_model1) #queries for model 1
  selected_queries2 <- select_queries(p_in_sample,p_gt_model2) #queries for model 2
  selected_queries3 <- select_queries(p_in_sample,p_gt_model3) #queries for model 3
  
  #out-of-sample error
  #benchmarks and model 1
  ag_model1 <- p_out_sample[which(rownames(p_out_sample)%in%selected_queries1),]
  mean_error1 <- apply(ag_model1,2,mean)
  median_error1 <- apply(ag_model1,2,median)
  mean_naiv12 <- mean_error1['naiv12']
  mean_arima <- mean_error1['arima']
  median_naiv12 <- median_error1['naiv12']
  median_arima <- median_error1['arima']
  
  
  #model 2
  ag_model2 <- p_out_sample[which(rownames(p_out_sample)%in%selected_queries2),]
  mean_error2 <- apply(ag_model2,2,mean)[p_gt_model2]
  median_error2 <- apply(ag_model2,2,median)[p_gt_model2]
  
  #model 3
  ag_model3 <- p_out_sample[which(rownames(p_out_sample)%in%selected_queries3),]
  mean_error3 <- apply(ag_model3,2,mean)[p_gt_model3]
  median_error3 <- apply(ag_model3,2,median)[p_gt_model3]
  
  mean_errors <- c(mean_naiv12,mean_arima,mean_error1[p_gt_model1],mean_error2,mean_error3)
  median_errors <- c(median_naiv12,median_arima,median_error1[p_gt_model1],median_error2,median_error3)
  output_errors <- t(data.frame(mean_errors,median_errors))
  #names(output_errors) <- c('naiv12','arima',p_gt_model1,p_gt_model2,p_gt_model3)
  
  #######DM results
  #coding for dm test data frame
  dm_coding <- function(p_model1,p_model2){
    return(paste(p_model1,'_vs_',p_model2,sep=''))
  }
  
  if(p_print_dm){
    print('DM results')
    #model 1
    first_better <- length(which(p_dm[,dm_coding(p_gt_model1,'arima')]<0.05))
    second_better <- length(which(p_dm[,dm_coding('arima',p_gt_model1)]<0.05))
    print(p_gt_model1)
    print(paste(length(selected_queries1),'series;',first_better,'times was',p_gt_model1,'superior;',second_better,'times was','arima','superior'))
    print('########################')
    
    #model 2
    first_better <- length(which(p_dm[,dm_coding(p_gt_model2,'arima')]<0.05))
    second_better <- length(which(p_dm[,dm_coding('arima',p_gt_model2)]<0.05))
    print(p_gt_model2)
    print(paste(length(selected_queries2),'series;',first_better,'times was',p_gt_model2,'superior;',second_better,'times was','arima','superior'))
    print('########################')
    
    #model 3
    first_better <- length(which(p_dm[,dm_coding(p_gt_model3,'arima')]<0.05))
    second_better <- length(which(p_dm[,dm_coding('arima',p_gt_model3)]<0.05))
    print(p_gt_model3)
    print(paste(length(selected_queries3),'series;',first_better,'times was',p_gt_model3,'superior;',second_better,'times was','arima','superior'))
    print('########################')
  }
  
  if(!p_print_dm){
    print('Mean and median errors')
    print(output_errors)
  }
  
  return()
  
}


###############analyse results
#inclusion decision according to part 4.2.1 of the thesis
include_germany <- which(inclusion_decision$gt_germany+inclusion_decision$gt_germany_modified1+inclusion_decision$gt_germany_modified2>1)
include_international <- which(inclusion_decision$gt_international+inclusion_decision$gt_international_modified1+inclusion_decision$gt_international_modified2>1)

######part 4.4.1 (comparison against arima)
#Worldwide GT data according to part 4.4.1 of the thesis
print('#####results part 4.4.1 worldwide GT data######')
#MAPE
print('MAPE')
results_region(in_sample_mape[include_international,],out_sample_mape[include_international,],dm_results[include_international,],'gt_international','gt_international_modified1','gt_international_modified2',F)
#MASE
print('MASE')
results_region(in_sample_mape[include_international,],out_sample_mase[include_international,],dm_results[include_international,],'gt_international','gt_international_modified1','gt_international_modified2',F)
#DM results
results_region(in_sample_mape[include_international,],out_sample_mape[include_international,],dm_results[include_international,],'gt_international','gt_international_modified1','gt_international_modified2',T)

#German GT data according to part 4.4.1 of the thesis
print('#####results part 4.4.1 German GT data######')
#MAPE
print('MAPE')
results_region(in_sample_mape[include_germany,],out_sample_mape[include_germany,],dm_results[include_germany,],'gt_germany','gt_germany_modified1','gt_germany_modified2',T)
#MASE
print('MASE')
results_region(in_sample_mape[include_germany,],out_sample_mase[include_germany,],dm_results[include_germany,],'gt_germany','gt_germany_modified1','gt_germany_modified2',F)
print('##########################################')

#comparison of different GT versions
print('#####results part 4.4.2 German GT data######')
print('MAPE')
compare_statistics(out_sample_mape[include_germany,],'gt_germany_modified2','gt_germany')
print('MASE')
compare_statistics(out_sample_mase[include_germany,],'gt_germany_modified2','gt_germany')

print('#####results part 4.4.2 worldwide GT data######')
print('MAPE')
compare_statistics(out_sample_mape[include_international,],'gt_international_modified2','gt_international')
print('MASE')
compare_statistics(out_sample_mase[include_international,],'gt_international_modified2','gt_international')

