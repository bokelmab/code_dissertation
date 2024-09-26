#wrapper functions
source('WrapperFunctions/time_calculations.R')

#An object to fit models and evaluate their performance
#By constructing the object, the benchmark models are fitted.
#gt models can be added using build_gt_model(). Performance measures of all fitted models
#can be provided using out_sample_mape(),,,

build_model <- function(p_arrival,p_end_tr=c(2015,8)){
  
  build_model_env <- environment()
  models <- NULL
  model_residuals <- NULL
  
  #########build ARIMA and seasonal naive benchmark models
  #divide arrival data in training and test sample
  ar_tr <- window(p_arrival,end=p_end_tr) #arrival training
  ar_ts <- window(p_arrival,start=increase_month(p_end_tr,1)) #arrival test
  
  #intervention variable for change of data collection in 2012
  inter12 <- ts(c(rep(100,36),rep(0,68)),start=start(p_arrival),frequency = 12)
  intervention <- window_as(inter12,ar_tr) #intervention in training set
  
  #ARIMA
  arima_order <- arimaorder(auto.arima(ar_tr,xreg=intervention,D=1)) #order selection
  arima_model <- Arima(ar_tr,xreg=intervention,order=arima_order[1:3],seasonal=arima_order[4:6],
                       include.drift=(arima_order[2]==0),method='CSS') #force inclusion of drift
  intervention <- inter12 #intervention on complete data set
  refit_arima <- Arima(p_arrival,xreg=intervention, model=arima_model)
  models[[2]] <- arima_model #save model
  names(models)[2] <- 'arima'
  model_residuals[[1]] <- refit_arima$residuals #save residuals
  names(model_residuals)[1] <- 'arima'
  
  #naiv12
  naiv12=ts(c(rep(0,12),p_arrival-stats::lag(p_arrival,-12)),start=start(p_arrival),
            frequency=frequency(p_arrival))
  model_residuals[[2]] <- naiv12 #save residuals
  names(model_residuals)[2] <- 'naiv12'
  
  output <- list(
    
    env = build_model_env,
    
    ##########build gt model
    build_gt_model = function(p_gt,p_lags,p_name){ 
      #p_name represents name of GT version
      #p_lags are the lags of GT data included into the model
      
      #divide gt data in training and test sample 
      gt_tr <- window(p_gt,end=p_end_tr)
      gt_ts <- window(p_gt,start=increase_month(p_end_tr,1))
      
      #construct training data with arrival data and (lagged) gt data
      data_tr <- lagged_predictors(ar_tr,gt_tr,p_lags)
      
      #fit gt model
      arrival <- data_tr[,1]
      input <- data.frame(data_tr[,-1],intervention=window_as(inter12,ar_tr)) #(lagged) GT data and intervention
      names(input)[-ncol(input)] <- names(data_tr)[-1]
      gt_arima_order <- arimaorder(auto.arima(arrival, xreg=input, D=1))
      gt_model <- Arima(arrival, xreg=input, order=gt_arima_order[1:3],
                        seasonal=gt_arima_order[4:6],include.drift=(gt_arima_order[2]==0),
                        method='CSS')
      
      #forecasting
      data_all <- lagged_predictors(p_arrival,p_gt,p_lags)
      arrival <- data_all[,1]
      input <- data.frame(data_all[,-1],intervention=inter12)
      names(input)[-ncol(input)] <- names(data_all)[-1]
      refit_gt <- Arima(arrival, xreg=input, model=gt_model)  
      
      #add gt_model to model list
      models <- get('models',build_model_env)
      models[[length(models)+1]] <- gt_model
      names(models)[length(models)] <- p_name
      assign('models',models,build_model_env)
      
      #add residuals to residual list
      model_residuals <- get('model_residuals',build_model_env)
      model_residuals[[length(model_residuals)+1]] <- refit_gt$residuals
      names(model_residuals)[length(model_residuals)] <- p_name
      assign('model_residuals',model_residuals,build_model_env)
      
    },
    
    get_models = function(p_model){
      models <- get('models',build_model_env)
      return(models[[p_model]])
    },
    
    get_residuals = function(){
      return(get('model_residuals',build_model_env))
    },
    
    #########rmse for all models in-sample 
    insample_rmse = function(){ 
      rmse <- function(p_residuals){
        #because some models do not fit the first years, we start at the fifth year
        p_residuals <- window(p_residuals,start=increase_month(start(arrival),48),end=p_end_tr)
        return(sqrt(mean(p_residuals^2)))
      }
      return(lapply(model_residuals,rmse))
    },
    
    #########rmse for all models out-of-sample
    outsample_rmse = function(){
      rmse <- function(p_residuals){
        p_residuals <- window(p_residuals,start=increase_month(p_end_tr,1))
        return(sqrt(mean(p_residuals^2)))
      }
      return(lapply(model_residuals,rmse))
    },
    
    #########mape for all models in-sample 
    insample_mape = function(){ 
      mape <- function(p_residuals){
        #because some models do not fit the first years, we start at the fifth year
        p_residuals <- window(p_residuals,start=increase_month(start(arrival),48),end=p_end_tr)
        p_arrival <- window(p_arrival,start=increase_month(start(p_arrival),48),end=p_end_tr)
        return(mean(abs(p_residuals/p_arrival)))
      }
      return(lapply(model_residuals,mape))
    },
    
    #########mape for all models out-of-sample
    outsample_mape = function(){
      mape <- function(p_residuals){
        p_residuals <- window(p_residuals,start=increase_month(p_end_tr,1))
        p_arrival <- window(p_arrival,start=increase_month(p_end_tr,1))
        return(mean(abs(p_residuals/p_arrival)))
      }
      return(lapply(model_residuals,mape))
    },
    
    #########mae for all models in-sample 
    insample_mae = function(){ 
      mae <- function(p_residuals){
        #because some models do not fit the first years, we start at the fifth year
        p_residuals <- window(p_residuals,start=increase_month(start(arrival),48),end=p_end_tr)
        return(mean(abs(p_residuals)))
      }
      return(lapply(model_residuals,mae))
    },
    
    #########mae for all models out-of-sample
    outsample_mae = function(){
      mae <- function(p_residuals){
        p_residuals <- window(p_residuals,start=increase_month(p_end_tr,1))
        return(mean(abs(p_residuals)))
      }
      return(lapply(model_residuals,mae))
    },
    
    #########mase for all models in-sample 
    insample_mase = function(){ 
      mase <- function(p_residuals){
        #because some models do not fit the first years, we start at the fifth year
        p_residuals <- window(p_residuals,start=increase_month(start(arrival),48),end=p_end_tr)
        nominator <- mean(abs(p_residuals))
        denominator <- mean(abs(window(model_residuals[['naiv12']],start=increase_month(start(arrival),48),end=p_end_tr)))
        return(nominator/denominator)
      }
      return(lapply(model_residuals,mase))
    },
    
    #########mase for all models out-of-sample
    outsample_mase = function(){
      mase <- function(p_residuals){
        p_residuals <- window(p_residuals,start=increase_month(p_end_tr,1))
        nominator <- mean(abs(p_residuals))
        denominator <- mean(abs(window(model_residuals[['naiv12']],start=increase_month(start(arrival),48),end=p_end_tr)))
        return(nominator/denominator)
      }
      return(lapply(model_residuals,mase))
    },
    
    #########residual diagnostics 
    #Diebold-Mariano Test (Test if model 1 performs better then model 2)
    diebold_mariano = function(p_model1,p_model2){
      
      error1 <- window(model_residuals[[p_model1]],start=increase_month(p_end_tr,1))
      error2 <- window(model_residuals[[p_model2]],start=increase_month(p_end_tr,1))
      mse <- dm.test(error1,error2,alternative='less',
                     power=2)
      mae <- dm.test(error1,error2,alternative='less',
                     power=1)
      return(list(mse=mse,mae=mae))
    }
  )
    
  
}