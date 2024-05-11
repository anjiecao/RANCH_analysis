


pooled_sd <- function(sd, n) {
  sqrt(sum((n - 1) * sd^2) / sum(n - 1))
}


connect_files <- function(trial_info, param_info, sim_data){
  
  linked_file <- trial_info %>% 
    distinct(trial_id, stim_id, fam_duration, violation_type) %>% 
    right_join(sim_data, by = c("trial_id", "stim_id")) %>% 
    select(param_id, index, fam_duration, violation_type, sample_n_mean, sample_n_std) %>% 
    right_join(param_info, by = c('param_id')) %>%
    mutate(trial_type = case_when(
      violation_type == "background" ~ "background", 
      violation_type == "identity" & (index == fam_duration ) ~ "deviant",
      TRUE ~ "background"))%>% 
    mutate(trial_number = index + 1) %>% 
    group_by(param_id, trial_type, trial_number) %>% 
    summarise(
      count = n(),
      #ssd = sum((20-1) * (sample_n_std ** 2)),
      #sd_sample = (ssd / 20 * count - count) ** .5, 
      sd_sample = pooled_sd(sample_n_std, 20),
      mean_sample = mean(sample_n_mean), 
      ub_sample = mean_sample + 1.96 * (sd_sample /((400 * count) ** .5)),
      lb_sample = mean_sample - 1.96 * (sd_sample /((400 * count) ** .5)),) %>% # missing SD here!
    ungroup() 
  
  return(linked_file)
  
}


connect_files_stim_type <- function(trial_info, param_info, sim_data){
  
  linked_file <- trial_info %>% 
    distinct(trial_id, stim_id, fam_duration, violation_type) %>% 
    right_join(sim_data, by = c("trial_id", "stim_id")) %>% 
    select(param_id, index, fam_duration, violation_type, sample_n_mean, sample_n_std) %>%
    mutate(trial_type = case_when(
      index == fam_duration ~ violation_type,
      TRUE ~ "fam" )) %>% 
    mutate(trial_number = index + 1) %>% 
    group_by(param_id, trial_type, trial_number) %>% 
    summarise(
      count = n(),
      #ssd = sum((20-1) * (sample_n_std ** 2)),
      #sd_sample = (ssd / 20 * count - count) ** .5, 
      sd_sample = pooled_sd(sample_n_std, 20),
      mean_sample = mean(sample_n_mean), 
      ub_sample = mean_sample + 1.96 * (sd_sample /((400 * count) ** .5)),
      lb_sample = mean_sample - 1.96 * (sd_sample /((400 * count) ** .5)),) %>% # missing SD here!
    ungroup() 
  
  return(linked_file)
  
}