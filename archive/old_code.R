# Subject-wise crossvalidation 
infant_data = infant_data %>% mutate(subj_id = paste0(experiment, subject_num))


# aligned EIG
infant_kfold_aligned_eig <- lapply(
  unique(aligned_data$param_id), 
  function(id){
    run_subject_wise_crossvalidation(p_id = id, 
                                     raw_d = infant_data, 
                                     sd = aligned_data)
  }
) %>% 
  bind_rows()


# resnet 50 EIG
infant_kfold_resnet50_eig <- lapply(
  unique(resnet50_data$param_id), 
  function(id){
    run_subject_wise_crossvalidation(p_id = id, 
                                     raw_d = infant_data, 
                                     sd = resnet50_data)
  }
) %>% 
  bind_rows()

# resnet 50 KL
infant_kfold_resnet50_KL <- lapply(
  unique(KL_data$param_id), 
  function(id){
    run_subject_wise_crossvalidation(p_id = id, 
                                     raw_d = infant_data, 
                                     sd = KL_data)
  }
) %>% 
  bind_rows()

# resnet 50 surprisal
infant_kfold_resnet50_surprisal_short <- lapply(
  unique(surprisal_data_short$param_id), 
  function(id){
    run_subject_wise_crossvalidation(p_id = id, 
                                     raw_d = infant_data, 
                                     sd = surprisal_data_short)
  }
) %>% 
  bind_rows()

infant_kfold_resnet50_surprisal_long <- lapply(
  unique(surprisal_data_long$param_id), 
  function(id){
    run_subject_wise_crossvalidation(p_id = id, 
                                     raw_d = infant_data, 
                                     sd = surprisal_data_long)
  }
) %>% 
  bind_rows()


# resnet 50 nonoise
infant_kfold_resnet50_nonoise <- lapply(
  unique(nonoise_data$param_id), 
  function(id){
    run_subject_wise_crossvalidation(p_id = id, 
                                     raw_d = infant_data, 
                                     sd = nonoise_data)
  }
) %>% 
  bind_rows()

# resnet 50 nolearning
infant_kfold_resnet50_nolearning <- lapply(
  unique(nolearning_data$param_id), 
  function(id){
    run_subject_wise_crossvalidation(p_id = id, 
                                     raw_d = infant_data, 
                                     sd = nolearning_data)
  }
) %>% 
  bind_rows()

infant_kfold_resnet50_linmodel <- run_crossvalidation_linmodel(infant_data)


saveRDS(infant_kfold_aligned_eig, here("cached_data/infant_kfold_aligned_eig.Rds"))
saveRDS(infant_kfold_resnet50_eig, here("cached_data/infant_kfold_resnet50_eig.Rds"))
saveRDS(infant_kfold_resnet50_KL, here("cached_data/infant_kfold_resnet50_KL.Rds"))
saveRDS(infant_kfold_resnet50_surprisal_short, here("cached_data/infant_kfold_resnet50_surprisal_short.Rds"))
saveRDS(infant_kfold_resnet50_surprisal_long, here("cached_data/infant_kfold_resnet50_surprisal_long.Rds"))

saveRDS(infant_kfold_resnet50_nonoise, here("cached_data/infant_kfold_resnet50_nonoise.Rds"))
saveRDS(infant_kfold_resnet50_nolearning, here("cached_data/infant_kfold_resnet50_nolearning.Rds"))

