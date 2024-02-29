
run_k_fold <- function(p_id, raw_d, sd){
  
  # get all the subjects 
  all_sbj <- unique(raw_d$prolific_id)
  n_folds = 10
  
  # currently testing on one set of parameter 
  rsquared <- vector(mode = "list", length = n_folds)
  rmse <- vector(mode = "list", length = n_folds)
  
  for (k in 1:n_folds){
    
    # separate the data in to train & test
    train_subject = sample(all_sbj, floor(length(all_sbj) * .9))
    train_data = summarize_behavioral_data(raw_d %>% filter(prolific_id %in% train_subject)) %>% ungroup() %>% left_join(sd %>% filter(param_id == p_id) %>% ungroup(), by = c("trial_number", "trial_type")) #%>% 
     # filter(trial_number != 11)
    test_data = summarize_behavioral_data(raw_d %>% filter(!prolific_id %in% train_subject)) %>% ungroup() %>% left_join(sd %>% filter(param_id == p_id) %>% ungroup(), by = c("trial_number", "trial_type")) #%>% 
      #filter(trial_number != 11)
    
    # fit the training set
    
    
    fitted_stats = colf_nlxb(mean_lt ~ mean_sample, data = train_data, lower = c(-Inf, 0.0000001))
    sample_slope = fitted_stats$coefficients["param_mean_sample"]
    sample_intercept = fitted_stats$coefficients["param_X.Intercept."]
    
    test_scaled <- test_data %>% 
      mutate(scaled_samples = mean_sample * sample_slope  + sample_intercept)
    
    # save results for each fold 
    rsquared[[k]] <- cor(test_scaled$mean_lt, test_scaled$scaled_samples)^2
    rmse[[k]] <- rmse(test_scaled$mean_lt, test_scaled$scaled_samples)
    
  }
  
  rsquared_average <- mean(unlist(rsquared))
  rmse_average <- mean(unlist(rmse))
  rsquared_sd = sd(unlist(rsquared))
  rmse_sd <- sd(unlist(rmse))
  
  
  return(
    tibble(
      "param_id" = p_id, 
      "mean_rsquared" = rsquared_average, 
      "mean_rmse" = rmse_average, 
      "ub_rsquared" = rsquared_average + 1.96 * (rsquared_sd / sqrt(n_folds)), 
      "lb_rsquared" = rsquared_average - 1.96 * (rsquared_sd / sqrt(n_folds)), 
      "ub_rmse" = rmse_average + 1.96 * (rmse_sd / sqrt(n_folds)), 
      "lb_rmse" = rmse_average - 1.96 * (rmse_sd / sqrt(n_folds))
      
    )
  )
  
}


link_dataset <- function(p_id, raw_d, sd){
  
 
  
  

  data = summarize_behavioral_data(raw_d) %>% ungroup() %>% left_join(sd %>% filter(param_id == p_id) %>% ungroup(), by = c("trial_number", "trial_type")) %>% 
    filter(trial_number != 11)
  
  
  return(
    data)
  
}

summarize_behavioral_data <- function(raw_d){
  
  bd_fit_summary <- raw_d %>% 
    mutate(fam_duration = exposure_duration) %>% 
    group_by(trial_number,fam_duration,trial_type) %>% 
    summarise(mean_lt = mean(total_rt), 
              count = n(),
              #ssd = sum((20-1) * (sample_n_std ** 2)),
              #sd_sample = (ssd / 20 * count - count) ** .5, 
              sd_lt = sd(total_rt),
              ub_lt = mean_lt + 1.96 * (sd_lt /(count ** .5)),
              lb_lt = mean_lt - 1.96 * (sd_lt /(count ** .5))) 
  
  return (bd_fit_summary)
}

summarize_behavioral_data_infants <- function(raw_d){
  
  bd_fit_summary <- raw_d %>% 
    group_by(fam_duration,test_type) %>% 
    summarise(mean_lt = mean(LT), 
              count = n(),
              #ssd = sum((20-1) * (sample_n_std ** 2)),
              #sd_sample = (ssd / 20 * count - count) ** .5, 
              sd_lt = sd(LT),
              ub_lt = mean_lt + 1.96 * (sd_lt /(count ** .5)),
              lb_lt = mean_lt - 1.96 * (sd_lt /(count ** .5))) 
  
  return (bd_fit_summary)
}

run_subject_wise_crossvalidation <- function(p_id, raw_d, sd){
  
  # get all the subjects 
  all_sbj <- unique(raw_d$subj_id)
  n_folds = length(all_sbj)
  
  # currently testing on one set of parameter 
  rsquared <- vector(mode = "list", length = n_folds)
  rmse <- vector(mode = "list", length = n_folds)
  
  for (k in 1:n_folds){
    
    # separate the data in to train & test
    train_subject = sample(all_sbj, n_folds - 1)
    train_data = summarize_behavioral_data_infants(raw_d %>% filter(subj_id %in% train_subject)) %>% ungroup() %>% left_join(sd %>% filter(param_id == p_id) %>% ungroup(), by = c("fam_duration", "test_type")) 
    test_data = summarize_behavioral_data_infants(raw_d %>% filter(!subj_id %in% train_subject)) %>% ungroup() %>% left_join(sd %>% filter(param_id == p_id) %>% ungroup(), by = c("fam_duration", "test_type")) 
    
    # fit the training set
    fitted_stats = colf_nlxb(mean_lt ~ mean_sample, data = train_data, lower = c(-Inf, 0.0000001))
    sample_slope = fitted_stats$coefficients["param_mean_sample"]
    sample_intercept = fitted_stats$coefficients["param_X.Intercept."]
    
    test_scaled <- test_data %>% 
      mutate(scaled_samples = mean_sample * sample_slope  + sample_intercept)
    
    # save results for each fold 
    rsquared[[k]] <- cor(test_scaled$mean_lt, test_scaled$scaled_samples)^2
    rmse[[k]] <- rmse(test_scaled$mean_lt, test_scaled$scaled_samples)
    
  }
  
  rsquared_average <- mean(unlist(rsquared))
  rmse_average <- mean(unlist(rmse))
  rsquared_sd = sd(unlist(rsquared))
  rmse_sd <- sd(unlist(rmse))
  
  
  return(
    tibble(
      "param_id" = p_id, 
      "mean_rsquared" = rsquared_average, 
      "mean_rmse" = rmse_average, 
      "ub_rsquared" = rsquared_average + 1.96 * (rsquared_sd / sqrt(n_folds)), 
      "lb_rsquared" = rsquared_average - 1.96 * (rsquared_sd / sqrt(n_folds)), 
      "ub_rmse" = rmse_average + 1.96 * (rmse_sd / sqrt(n_folds)), 
      "lb_rmse" = rmse_average - 1.96 * (rmse_sd / sqrt(n_folds))
      
    )
  )
  
}
