run_joint_fit <- function(adult_d, infant_d, adult_sim, infant_sim){
  
  adult_data_to_fit <- summarize_behavioral_data(adult_bd) %>% 
    ungroup() %>% left_join(adult_resnet50_eig, by = c("trial_number", "trial_type")) 
  
  infant_data_to_fit <- summarize_behavioral_data_infants(infant_d) %>% 
    ungroup() %>% 
    left_join(infant_resnet50_eig %>% ungroup(), by = c("fam_duration", "test_type")) 
  
  joint_scale_dataset <- adult_data_to_fit %>% 
    mutate(type = "adults") %>% 
    #select(type, mean_lt, mean_sample) %>% 
    mutate(mean_lt = mean_lt / 1000) %>% 
    bind_rows(
      infant_data_to_fit %>% 
        mutate(type = "infants")  
      #select(type, mean_lt, mean_sample)
    ) 
  
  fit_df = tibble()
  
  param_infants = unique(infant_sim$param_id)
  param_adults = unique(adult_sim$param_id)
  
  counter = 0
  
  fit_list = vector("list", length = length(param_infants) * length(param_adults))
  
  for (p_i in param_infants) {
    
    print("infant_param")
    print(p_i)
    
    for (p_a in param_adults) {
      
      counter = counter + 1
      
      # subset data
      cur_data = joint_scale_dataset %>% 
        filter((param_id == p_i & type == 'infants') | 
                 (param_id == p_a & type == 'adults'))
      
      fitted_stats = lm(mean_lt ~ 0 + mean_sample, data = cur_data)
      sample_slope = fitted_stats$coefficients["mean_sample"]
      sample_intercept = 0 #fitted_stats$coefficients["(Intercept)"]
      
      scaled_data <- cur_data %>% 
        mutate(scaled_samples = mean_sample * sample_slope  + sample_intercept)
      
      #calculate the r2 and rmse
      rsquared <- cor(scaled_data$mean_lt, scaled_data$scaled_samples)^2
      rmse <- rmse(scaled_data$mean_lt, scaled_data$scaled_samples)
      
      fit_list[[counter]] = list(p_i = p_i, p_a = p_a, rsquared = rsquared, rmse = rmse)
      
    }
  }
  
  
  fit_df = bind_rows(fit_list)
  
  return(fit_df)
}