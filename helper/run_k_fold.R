
run_all_kfold <- function(sim_data, behavioral_data){
  
  full_k_fold_res <- lapply(
    unique(sim_data$param_id), 
    function(id){
      run_k_fold(p_id = id, 
                 raw_d = behavioral_data, 
                 sd = sim_data)
    }
  ) %>% 
    bind_rows()
  
  return (full_k_fold_res)
  
}

split_vector_into_chunks <- function(vec){
  n <- length(vec) # Includes NA in the count
  base_size <- n %/% 10
  extra_elements <- n %% 10
  
  sizes <- rep(base_size, 10)
  sizes[1:extra_elements] = sizes[1:extra_elements] + 1
  
  chunks <- vector("list", length = 10)
  start_index <- 1
  for (i in 1:10) {
    end_index <- start_index + sizes[i] - 1
    chunks[[i]] <- vec[start_index:end_index]
    start_index <- end_index + 1
  }
  
  return(chunks)
}

run_all_full_data_fit <- function(sim_data, behavioral_data){
  
  full_data_fit_res <- lapply(
    unique(sim_data$param_id), 
    function(id){
      run_full_data_fit(p_id = id, 
                 raw_d = behavioral_data, 
                 sd = sim_data)
    }
  ) %>% 
    bind_rows()
  
  return (full_data_fit_res)
  
}


run_full_data_fit <- function(p_id, raw_d, sd){
  
  data_to_fit <- summarize_behavioral_data(raw_d) %>% 
    ungroup() %>% 
    left_join(sd %>% filter(param_id == p_id) %>% ungroup(), by = c("trial_number", "trial_type")) 
  
  fitted_stats = colf_nlxb(mean_lt ~ mean_sample, data = data_to_fit, lower = c(-Inf, 0.0000001))
  sample_slope = fitted_stats$coefficients["param_mean_sample"]
  sample_intercept = fitted_stats$coefficients["param_X.Intercept."]
  
  scaled_data <- data_to_fit %>% 
    mutate(scaled_samples = mean_sample * sample_slope  + sample_intercept)
  
  rsquared <- cor(scaled_data$mean_lt, scaled_data$scaled_samples)^2
  rmse <- rmse(scaled_data$mean_lt, scaled_data$scaled_samples)
 
  
  return(
    tibble(
      "param_id" = p_id, 
      "mean_rsquared" = rsquared, 
      "mean_rmse" = rmse
    )
  )
  
}

run_k_fold <- function(p_id, raw_d, sd){
  
  # get all the subjects 
  all_sbj <- unique(raw_d$prolific_id)
  all_sbj <- sample(all_sbj)
  n_folds = 10
  
  
  folds = split_vector_into_chunks(all_sbj)

  
  # currently testing on one set of parameter 
  rsquared <- vector(mode = "list", length = n_folds)
  rmse <- vector(mode = "list", length = n_folds)
  
  for (k in 1:n_folds){
    
    # separate the data in to train & test
    train_subject = folds[[k]]
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
    
    # generate index
    ind = rep(TRUE, n_folds)
    ind[k] = FALSE
    
    # separate the data in to train & test
    train_subject = all_sbj[ind]
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

run_splithalf_crossvalidation <- function(p_id, raw_d, sd){
  
  # get all the subjects 
  all_sbj <- unique(raw_d$subj_id)
  n_folds = 2
  
  # currently testing on one set of parameter 
  rsquared <- vector(mode = "list", length = n_folds)
  rmse <- vector(mode = "list", length = n_folds)
  
  # random indices
  ind <- sample(c(TRUE, FALSE), length(all_sbj), replace=TRUE, prob=c(0.5, 0.5))

  # create two equal sets
  set1 = all_sbj[ind]
  set2 = all_sbj[!ind]
  
  for (k in 1:n_folds){
    
    # separate the data in to train & test
    if (k==1) {
      train_subject = set1
    }
    else {
      train_subject = set2
    }
    
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

run_splithalf_linmodel <- function(raw_d) {
  # get all the subjects 
  all_sbj <- unique(raw_d$subj_id)
  n_folds = 2
  
  # currently testing on one set of parameter 
  rsquared <- vector(mode = "list", length = n_folds)
  rmse <- vector(mode = "list", length = n_folds)
  
  # random indices
  ind <- sample(c(TRUE, FALSE), length(all_sbj), replace=TRUE, prob=c(0.5, 0.5))
  
  # create two equal sets
  set1 = all_sbj[ind]
  set2 = all_sbj[!ind]
  
  for (k in 1:n_folds){
    
    # separate the data in to train & test
    if (k==1) {
      train_subject = set1
    }
    else {
      train_subject = set2
    }
    
    train_data = summarize_behavioral_data_infants(raw_d %>% filter(subj_id %in% train_subject)) 
    test_data = summarize_behavioral_data_infants(raw_d %>% filter(!subj_id %in% train_subject))  
    
    
    # fit linear model
    linear_model = tidy(lm(mean_lt ~ fam_duration * test_type - test_type, data = train_data))
    
    intercept_linmodel = linear_model$estimate[1]
    fam_duration_linmodel = linear_model$estimate[2]
    interaction_linmodel = linear_model$estimate[3]
    
    test_scaled = test_data %>%
      mutate(samples_linear_model = intercept_linmodel + fam_duration_linmodel * fam_duration + interaction_linmodel * ifelse(test_type == 'nov', 1, 0) * fam_duration)
    
    # save results for each fold 
    rsquared[[k]] <- cor(test_scaled$mean_lt, test_scaled$samples_linear_model)^2
    rmse[[k]] <- rmse(test_scaled$mean_lt, test_scaled$samples_linear_model)
    
  }
  
  rsquared_average <- mean(unlist(rsquared))
  rmse_average <- mean(unlist(rmse))
  rsquared_sd = sd(unlist(rsquared))
  rmse_sd <- sd(unlist(rmse))
  
  
  return(
    tibble(
      "mean_rsquared" = rsquared_average, 
      "mean_rmse" = rmse_average, 
      "ub_rsquared" = rsquared_average + 1.96 * (rsquared_sd / sqrt(n_folds)), 
      "lb_rsquared" = rsquared_average - 1.96 * (rsquared_sd / sqrt(n_folds)), 
      "ub_rmse" = rmse_average + 1.96 * (rmse_sd / sqrt(n_folds)), 
      "lb_rmse" = rmse_average - 1.96 * (rmse_sd / sqrt(n_folds))
      
    )
  )
  }
  
run_crossvalidation_linmodel <- function(raw_d) {
  
  # get all the subjects 
  all_sbj <- unique(raw_d$subj_id)
  n_folds = length(all_sbj)
  
  # currently testing on one set of parameter 
  rsquared <- vector(mode = "list", length = n_folds)
  rmse <- vector(mode = "list", length = n_folds)
  
  for (k in 1:n_folds){
    
    # generate index
    ind = rep(TRUE, n_folds)
    ind[k] = FALSE
    
    # separate the data in to train & test
    train_subject = all_sbj[ind]
    train_data = summarize_behavioral_data_infants(raw_d %>% filter(subj_id %in% train_subject)) 
    test_data = summarize_behavioral_data_infants(raw_d %>% filter(!subj_id %in% train_subject))  
    
    # fit linear model
    linear_model = tidy(lm(mean_lt ~ fam_duration * test_type - test_type, data = train_data))
    
    intercept_linmodel = linear_model$estimate[1]
    fam_duration_linmodel = linear_model$estimate[2]
    interaction_linmodel = linear_model$estimate[3]
    
    test_scaled = test_data %>%
      mutate(samples_linear_model = intercept_linmodel + fam_duration_linmodel * fam_duration + interaction_linmodel * ifelse(test_type == 'nov', 1, 0) * fam_duration)
    
    # save results for each fold 
    rsquared[[k]] <- cor(test_scaled$mean_lt, test_scaled$samples_linear_model)^2
    rmse[[k]] <- rmse(test_scaled$mean_lt, test_scaled$samples_linear_model)
    
  }
  
  rsquared_average <- mean(unlist(rsquared))
  rmse_average <- mean(unlist(rmse))
  rsquared_sd = sd(unlist(rsquared))
  rmse_sd <- sd(unlist(rmse))
  
  
  return(
    tibble(
      "mean_rsquared" = rsquared_average, 
      "mean_rmse" = rmse_average, 
      "ub_rsquared" = rsquared_average + 1.96 * (rsquared_sd / sqrt(n_folds)), 
      "lb_rsquared" = rsquared_average - 1.96 * (rsquared_sd / sqrt(n_folds)), 
      "ub_rmse" = rmse_average + 1.96 * (rmse_sd / sqrt(n_folds)), 
      "lb_rmse" = rmse_average - 1.96 * (rmse_sd / sqrt(n_folds))
      
    )
  )
  
}
  
