viz_top_param <- function(sim_data, infant_data, top_param) {
  
  combined_df = sim_data %>% mutate(fam_duration = as.integer(trial_number) - 1) %>%
    left_join(grouped_infant_means, by = c("test_type", "fam_duration")) %>%
    filter(!is.na(LT)) 
  
  grouped_stats <- combined_df %>%
    group_by(param_id) %>%
    do(tidy(lm(LT ~ mean_sample, data = .))) %>%
    select(-c(std.error, statistic, p.value)) %>%
    ungroup() %>%
    spread(term, estimate)
  
  # merge with intercept and slope and scale the samples
  combined_df <- combined_df %>%
    left_join(grouped_stats, by = c("param_id")) %>% 
    dplyr::rename(intercept = "(Intercept)", slope = mean_sample.y, n_samples = mean_sample.x) %>% 
    mutate(scaled_samples = n_samples * slope + intercept,
           scaled_ub_sample = ub_sample * slope  + intercept, 
           scaled_lb_sample = lb_sample * slope + intercept) %>% select(-n_samples)
  
  top_sim_result = combined_df %>% filter(param_id == top_param) 
  
  plot_df = top_sim_result %>% dplyr::rename(`Infant behavior` = LT, `RANCH (infants)` = scaled_samples) %>%
    pivot_longer(cols = c("RANCH (infants)", "Infant behavior"), names_to = "value_type", values_to = "value") %>% 
    mutate(value_type = factor(value_type, levels = c("Infant behavior", "RANCH (infants)")),
           test_type = ifelse(test_type == 'nov', 'Novel', 'Familiar')) 
  
  plot = ggplot(plot_df, aes(x = fam_duration, y = value, color = test_type, group = test_type)) + geom_smooth(method = "lm", size = 2) +  geom_point(size = 3)  + ylab('Looking time (s)') + xlab('fam duration') +  scale_x_continuous(breaks = c(0,1,2,3,4,6,8,9)) + xlab('Exposure duration') + labs(color = 'Test type') + theme(legend.position = "bottom", legend.text = element_text(size = 25)) + facet_grid(~value_type) + theme_classic(30) + scale_color_manual(values = c("blue", "red"))
  
  print(plot)
  
  return(plot_df %>% filter(value_type == 'RANCH (infants)'))
}

renaming <- function(df) {
  df <- df %>%
    group_by(param_id, trial_type) %>%
    filter(all(2:10 %in% trial_number)) %>% # only keep parameters for which all fam durations have been computed
    ungroup() %>%
    rename(test_type = trial_type) %>%
    mutate(fam_duration = trial_number - 1,
           test_type = case_when(test_type == 'background' ~ 'fam', test_type == 'deviant' ~ 'nov')) 
  return(df)
}

get_param_stats <- function(df){
  df = df %>% filter(!is.na(mean_rsquared)) %>% group_by(param_id) %>% 
    summarize(r2 = mean(mean_rsquared), rmse = mean(mean_rmse), 
              r2_lb = quantile(mean_rsquared, 0.025), r2_ub = quantile(mean_rsquared, 0.975),
              rmse_lb = quantile(mean_rmse, 0.025), rmse_ub = quantile(mean_rmse, 0.975))
  return(df)
}



scale_all_adult_sim <- function(sim_data, raw_d){
  
  full_data_scaled <- lapply(
    unique(sim_data$param_id), 
    function(id){
      scale_adult_sim(p_id = id, 
                        sim_data = sim_data, 
                        raw_d = raw_d) %>% 
        mutate(param_id = id)
    }
  ) %>% 
    bind_rows()
  
  return (full_data_scaled)
  
}

scale_all_infant_sim <- function(sim_data, raw_d){
  
  full_data_scaled <- lapply(
    #seq(1, 53, 1),
    seq(1, 161, 1),  #hardcoded for now since it only runs the first 162 params 
    #unique(sim_data$param_id), 
    function(id){
      print(id)
      scale_infant_sim(p_id = id, 
                      sim_data = sim_data, 
                      grouped_infant_means = grouped_infant_means) %>% 
        mutate(param_id = id)
    }
  ) %>% 
    bind_rows()
  
  return (full_data_scaled)
  
}

scale_adult_sim <- function(p_id, sim_data, raw_d){
  
  data_to_fit <- summarize_behavioral_data(raw_d) %>% 
    ungroup() %>% 
    left_join(sim_data %>% filter(param_id == p_id) %>% ungroup(), by = c("trial_number", "trial_type")) 
  
  fitted_stats = colf_nlxb(mean_lt ~ mean_sample, data = data_to_fit, lower = c(-Inf, 0.0000001))
  sample_slope = fitted_stats$coefficients["param_mean_sample"]
  sample_intercept = fitted_stats$coefficients["param_X.Intercept."]
  
  scaled_data <- data_to_fit %>% 
    mutate(scaled_samples = mean_sample * sample_slope  + sample_intercept, 
           scaled_lb = lb_sample * sample_slope + sample_intercept, 
           scaled_ub = ub_sample * sample_slope + sample_intercept)
  
  return(scaled_data)

  
}


scale_infant_sim <- function(p_id, sim_data, grouped_infant_means) {
  
  
  combined_df = sim_data %>% filter(param_id == p_id) %>% 
    mutate(fam_duration = as.integer(trial_number) - 1) %>%
    rename(test_type = trial_type) %>% 
    mutate(test_type = if_else(test_type == "background", "fam", "nov")) %>% 
    left_join(grouped_infant_means, by = c("test_type", "fam_duration")) %>%
    filter(!is.na(LT)) 
  
  stats <- combined_df %>%
    do(broom::tidy(lm(LT ~ mean_sample, data = .))) %>%
    select(-c(std.error, statistic, p.value)) %>%
    ungroup() %>%
    spread(term, estimate) %>% 
    rename(intercept = "(Intercept)", slope = mean_sample)
  
  # merge with intercept and slope and scale the samples
  combined_df <- combined_df %>%
    rename(n_samples = mean_sample) %>% 
    mutate(scaled_samples = n_samples * stats$slope + stats$intercept,
           scaled_ub_sample = ub_sample * stats$slope  + stats$intercept, 
           scaled_lb_sample = lb_sample * stats$slope + stats$intercept) %>% select(-n_samples)
  
  

  return(combined_df)
}



