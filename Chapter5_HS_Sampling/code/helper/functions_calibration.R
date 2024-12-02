create_Platt_outcome_prediction <- function(p_data, p_rf_model) {
  
  rf_model <- p_rf_model
  p_data$pred_prob <- predict(rf_model, data = p_data)$predictions
  glm_cal <- glm(Y ~ pred_prob, data = p_data, family = binomial)
  
  Platt_outcome_function <- function(p_data) {
    data_new <- data.frame(pred_prob = predict(rf_model, data = p_data)$predictions)
    return(predict(glm_cal, newdata = data_new, type = "response"))
  }
  
  return(Platt_outcome_function)  
}

create_isotonic_outcome_prediction <- function(p_data, p_rf_model){
  
  rf_model <- p_rf_model
  p_data$pred_prob <- predict(rf_model, data = p_data)$predictions
  isoreg_model <- isoreg(p_data$pred_prob, p_data$Y)
  isoreg_func <- as.stepfun(isoreg_model)
  
  istonic_outcome_function <- function(p_data) {
    pred_prob = predict(rf_model, data = p_data)$predictions
    return(isoreg_func(pred_prob))
  }
}