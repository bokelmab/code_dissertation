#required files
source('Forecasting/build_model.R') #forecasting and forecast evaluation
source('WrapperFunctions/cross_correlations.R') #pre-whitened cross-correlation
source('Forecasting/gt_analysis.R') #inclusion decision and order selection

env_workflow_competition <- environment()

gt_regions <- c('germany','international') 
versions <- c('','_modified1','_modified2') #original, mean divided and trend divided
gt_data_path <- 'Data/GT_Data/'

#load GT and arrival data
#Each version of GT data has its own data frame
assign('gt_germany_data',readRDS(paste(gt_data_path,'germany',"/gt_data_germany_final.RDS",sep='')),env_workflow_competition)
assign('gt_international_data',readRDS(paste(gt_data_path,'international',"/gt_data_international_final.RDS",sep='')),env_workflow_competition)
assign('arrival_data',readRDS("Data/Arrival_Data/arrival_data_matched.RDS"),env_workflow_competition)
for(i in 1:length(gt_regions)){
  for(j in 2:length(versions)){
    path <- paste(gt_data_path,'modified','/',gt_regions[i],versions[j],'.RDS',sep='')
    name <- paste('gt_',gt_regions[i],'_data',versions[j],sep='')
    assign(name,readRDS(path),env_workflow_competition)
  }
}


##########data frames to save the results

#query names (of included time series)
names_included <- c()

#error metrics
in_sample_rmse <- NULL
out_sample_rmse <- NULL
in_sample_mape <- NULL
out_sample_mape <- NULL
in_sample_mase <- NULL
out_sample_mase <- NULL

#inclusion decision
inclusion_decision <- NULL

#cross-correlation at all lags
for(i in 1:length(gt_regions)){
  for(j in 1:length(versions)){
    assign(paste('cross_correlations_',gt_regions[i],versions[j],sep=''),
           NULL,env_workflow_competition)
  }
}

#gt coefficient information
for(i in 1:length(gt_regions)){
  for(j in 1:length(versions)){
    assign(paste('coefficient_information_',gt_regions[i],versions[j],sep=''),
           NULL,env_workflow_competition)
  }
}

#diebold-mariano results
dm_results <- NULL

#exclusion of time gt time series where the respective holiday regions tourist data does not
#meet the specified requirements (No missing values, no change in name)
columns_to_delete <- which(is.na(arrival_data[1,]))

arrival_data <- arrival_data[,-columns_to_delete]
for(i in 1:length(gt_regions)){
  for(j in 1:length(versions)){
    name <- paste('gt_',gt_regions[i],'_data',versions[j],sep='')
    data <- get(name,env_workflow_competition)
    data <- data[,-columns_to_delete]
    assign(name,data,env_workflow_competition)
  }
}


##########Start of forecasting competition 
#The following for loop estimates for each pair of a holiday region and corresponding GT
#data a forecasting model. According to parts 4.2.1, 4.4.1 and 4.4.2 of the thesis some
#selection needs to be done to get the results of this thesis.

#PLEASE NOTE that the for loop takes arround 3 hours. Since the results of the for loop are safed
#in Forecasting/Results it can be skipped and the part "Analyse competition results" can be
#run next.

#The columns of each of the GT data frames (with the different versions of GT data) correspond
#to the 269 collected queries. Column i of the arrival_data frame includes the arrival data
#of the holiday region for which query i is a candidate query
#For each of the 269 queries the pre-whitened cross-correlation at lag 0 is calculated and
#forecasting models are build and evaluated.
#The file analyse_competition_results.R applies the query selection criterion described in
#part 4.2.1 of the thesis and sumerrizes the results.
for(i in 1:ncol(gt_germany_data)){ 
                                   
  ##########load data
  #get arrival data
  name_query <- names(gt_germany_data)[i]
  arrival <- arrival_data[,i]
  
  #get all versions of GT data for column i 
  for(j in 1:length(gt_regions)){
    for(k in 1:length(versions)){
      name <- paste('gt_',gt_regions[j],versions[k],sep='')
      data <- get(paste('gt_',gt_regions[j],'_data',versions[k],sep=''),env_workflow_competition)
      assign(name,data[,i],env_workflow_competition)
    }
  }
  
  ##########1. step: Inclusion decision and Gt lag order selection
  for(j in 1:length(gt_regions)){
    for(k in 1:length(versions)){
      name <- paste('gt_analysis_',gt_regions[j],versions[k],sep='')
      name_ts <- paste('gt_',gt_regions[j],versions[k],sep='')
      assign(name,gt_analysis(get(name_ts,env_workflow_competition),arrival))
    }
  }
  
  #forecast only if there is GT data with significant (pre-whitened) cross-correlation 
  individual_decision <- c()
  for(j in 1:length(gt_regions)){
    for(k in 1:length(versions)){
      data <- get(paste('gt_analysis_',gt_regions[j],versions[k],sep=''),env_workflow_competition)
      #significant positive pcc at lag 0 for specific GT version?
      individual_decision <- c(individual_decision,data[['include']]) 
    }
  }
  
  include <- (sum(individual_decision)>0) #at least one version with significant, positive pcc?
  #This is a pre-selection to save runtime. In the file analyse_competition_results the strict
  #selection criterion described in part 4.1.2 of the thesis is applied
  
  ##########2. step: Model building and forecasting
  if(include){
    model <- build_model(arrival) #benchmark models
    
    #build GT models for all GT versions
    for(j in 1:length(gt_regions)){ 
      for(k in 1:length(versions)){
        data <- get(paste('gt_analysis_',gt_regions[j],versions[k],sep=''),env_workflow_competition)
        name <- paste('gt_',gt_regions[j],versions[k],sep='')
        gt_ts <- get(name,env_workflow_competition)
        model$build_gt_model(gt_ts,data[['lags']],name) #GT model for specific GT version
      }
    }
    
  }
  
  #########3. step: Get competition results for query i
  if(include){
    names_included <- c(names_included,names(gt_germany_data)[i]) #name of query
    
    #error metrics
    in_sample_rmse <- rbind(in_sample_rmse,unlist(model$insample_rmse()))
    out_sample_rmse <- rbind(out_sample_rmse,unlist(model$outsample_rmse()))
    in_sample_mape <- rbind(in_sample_mape,unlist(model$insample_mape()))
    out_sample_mape <- rbind(out_sample_mape,unlist(model$outsample_mape()))
    in_sample_mase <- rbind(in_sample_mase,unlist(model$insample_mase()))
    out_sample_mase <- rbind(out_sample_mase,unlist(model$outsample_mase()))
    
    #inclusion decision (of individual gt data type)
    inclusion_decision <- rbind(inclusion_decision,individual_decision)
    
    #(pre-whitened) cross-correlations
    for(j in 1:length(gt_regions)){
      for(k in 1:length(versions)){
        name <- paste('cross_correlations_',gt_regions[j],versions[k],sep='')
        name_new <- paste('gt_analysis_',gt_regions[j],versions[k],sep='')
        data <- get(name,env_workflow_competition)
        data_new <- get(name_new,env_workflow_competition)
        data <- rbind(data,data_new[['cross_correlations']])
        assign(name,data)
      }
    }
    
    #coefficient information
    for(j in 1:length(gt_regions)){
      for(k in 1:length(versions)){
        name <- paste('coefficient_information_',gt_regions[j],versions[k],sep='')
        data <- get(name,env_workflow_competition)
        coef <- model$get_models(paste('gt_',gt_regions[j],versions[k],sep=''))$coef['gt']
        sd_coef <- model$get_models(paste('gt_',gt_regions[j],versions[k],sep=''))$var.coef['gt','gt']
        data <- rbind(data,c(coef,sqrt(sd_coef)))
        assign(name,data,env_workflow_competition)
      }
    }
    
    #diebold-mariano results
    names_models <- names(model$get_residuals())
    new_dm_results <- c()
    for(j in (1:length(names_models))){
      for(k in (1:length(names_models))[-j]){
        new_dm_results <- c(new_dm_results,model$diebold_mariano(names_models[j],names_models[k])$mae$p.value)
      }
    }
    dm_results <- rbind(dm_results,new_dm_results)
                                          
  }
  
  #provide intermediate information
  print(i)
  if(i%%10==0){
    print('MAPE')
    print(apply(out_sample_mape,2,mean))
    print('MASE')
    print(apply(out_sample_mase,2,mean))
    print(paste('Number of included series: ',nrow(out_sample_mape)))
    print('The results are just an oriantiation. Later query and model selection according to parts 4.2.1, 4.4.1 and 4.4.2 of the thesis has to be done.')
  }
  
}

#name data frames with competition information and save results
add_rownames <- function(p_frame_names,p_rownames){
  for(i in 1:length(p_frame_names)){
    frame <- get(p_frame_names[i],env_workflow_competition)
    rownames(frame) <- NULL
    frame <- data.frame(frame)
    rownames(frame) <- p_rownames
    assign(p_frame_names[i],frame,env_workflow_competition)
  }
}

add_colnames <- function(p_frame_names,p_colnames){
  for(i in 1:length(p_frame_names)){
    frame <- get(p_frame_names[i],env_workflow_competition)
    frame <- data.frame(frame)
    names(frame) <- p_colnames
    assign(p_frame_names[i],frame,env_workflow_competition)
  }
}

list_gt_types <- c()
for(i in 1:length(gt_regions)){
  for(j in 1:length(versions)){
    list_gt_types <- c(list_gt_types,paste(gt_regions[i],versions[j],sep=''))
  }
}

add_rownames(c('in_sample_rmse','out_sample_rmse','in_sample_mape','out_sample_mape',
               'in_sample_mase','out_sample_mase','inclusion_decision',
               paste('cross_correlations_',list_gt_types,sep=''),
               paste('coefficient_information_',list_gt_types,sep='')),names_included)


add_colnames(paste('coefficient_information_',list_gt_types,sep=''),c('gt_coef','sd'))
add_colnames(paste('cross_correlations_',list_gt_types,sep=''),as.character((-17):17))
add_colnames('inclusion_decision',paste('gt_',list_gt_types,sep=''))

#dm results
dm_names <- c()
for(j in (1:length(names_models))){
  for(k in (1:length(names_models))[-j]){
    dm_names <- c(dm_names,paste(names_models[j],'_vs_',names_models[k],sep=''))
  }
}
row.names(dm_results) <- NULL
add_colnames('dm_results',dm_names)

###########save results
saveRDS(in_sample_mape,'Forecasting/Results/Errors/in_sample_mape.RDS')
saveRDS(in_sample_mase,'Forecasting/Results/Errors/in_sample_mase.RDS')
saveRDS(in_sample_rmse,'Forecasting/Results/Errors/in_sample_rmse.RDS')
saveRDS(out_sample_mape,'Forecasting/Results/Errors/out_sample_mape.RDS')
saveRDS(out_sample_mase,'Forecasting/Results/Errors/out_sample_mase.RDS')
saveRDS(out_sample_rmse,'Forecasting/Results/Errors/out_sample_rmse.RDS')
saveRDS(inclusion_decision,'Forecasting/Results/Errors/inclusion_decision.RDS')

for(i in 1:length(list_gt_types)){
  #cross correlation
  name <- paste('cross_correlations_',list_gt_types[i],sep='')
  data <- get(name,env_workflow_competition)
  saveRDS(data,paste('Forecasting/Results/Cross_Correlations/',name,'.RDS',sep=''))
  
  #coefficient information
  name <- paste('coefficient_information_',list_gt_types[i],sep='')
  data <- get(name,env_workflow_competition)
  saveRDS(data,paste('Forecasting/Results/Coefficient_Information/',name,'.RDS',sep=''))
}

saveRDS(dm_results,'Forecasting/Results/DM_Results/dm_results.RDS')

#############Analyse competition results
source('Forecasting/analyse_competition_results.R')

#about the warnings
print('The warnings about NaNs in sqrt(diag(x$var.coef)) are due to numerical problems of calculating standard errors for a DLM which is close to being non-stationary. The estimated model is however valid (see http://r.789695.n4.nabble.com/warnings-in-ARMA-with-other-regressor-variables-td4666039.html) ')
