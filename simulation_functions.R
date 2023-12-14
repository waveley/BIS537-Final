# ####################################
# 
# Program: simulation_functions.R
#
# Description: Contains functions to run full simulation
#
# Returns: List of dataframes according to number run
#          with covariates, propensity scores, potential outcomes,
#          observed outcome
# 
# ####################################

library(survival)
library(tidyverse)
library(riskRegression)
library(tidyverse)
library(prodlim)
library(progress)
library(splitstackshape)
library(beepr)

source("misc_functions.R")
source("propensity_score_functions.R")
source("outcome_functions.R")

gen_data <- function(n, n_df, betas, alphas, gammas){
  cov_list <- gen_cur_covs(n = n*n_df)
  cov_tib <- 
    tibble(
      x1 = cov_list[[1]], 
      x2 = cov_list[[2]], 
      x3 = cov_list[[3]]
      )
  
  cov_tib <- gen_prop_scores(betas, covs = cov_tib)
  out_tib <- gen_outcomes(cov_tib, a = alphas)
  out_cens_time <- gen_censor(out_tib, g = gammas) %>%
    mutate(
      obs_time = ifelse(surv_time < cens_time, surv_time, cens_time)
    )
  
  out_df_list <- list()
  for(i in 1:n_df){
    out_df_list[[i]] <- out_cens_time %>% slice(((i-1)*n + 1):(i*n))
  }
  
  return(out_df_list)
}

 gen_res <- function(df){
   ### STEP 0: UNADJUSTED ESTIMATES ###
   # unadjusted_km_est <- 
   #    survfit(
   #      Surv(obs_time, 1-cens_ind) ~ trt, 
   #      data = df)
   #  
   #  cur_km_est_unadj <- 
   #    unadjusted_km_est %>% 
   #    summary(times = 5)
    
    # unadjusted_cox_est <-
    #   coxph(
    #     Surv(obs_time, 1 - cens_ind) ~ strata(trt) + x1 + x2 + x3,
    #     data = df
    #   )
    # 
    # cur_cox_est_unadj <- 
    #   survfit(unadjusted_cox_est) %>% 
    #   summary(times = 5)
    
    #cur_cox_est_unadj
    #cur_km_est_unadj
    
    ### STEP 1: IPTW ###
    iptw_reg <- 
      glm(
        trt ~ x1 + x2 + x3, 
        family = "quasibinomial", 
        data = df)
    
    df$prop_score <- iptw_reg$fitted.values
    
    df <- 
      df %>% mutate(
        trt_weight = 
          ifelse(trt == 1, 
                 1/prop_score, 
                 1/(1 - prop_score))
      )
    
    df1 <- df
  
    ### STEP 3: Treatment-weighted Estimates ###
    
    trt_adjusted_km_est <- 
      survfit(
        Surv(obs_time, 1-cens_ind) ~ trt, 
        data = df,
        weights = trt_weight)
    
    
    cur_km_est_trt_adj <- 
      trt_adjusted_km_est %>% 
      summary(times = 5)
     
    # 
    # trt_adjusted_cox_est <-
    #   coxph(
    #     Surv(obs_time, 1 - cens_ind) ~ strata(trt) + x1 + x2 + x3,
    #     data = df,
    #     weights = trt_weight
    #   )
    # 
    # cur_cox_est_trt_adj <- 
    #   survfit(trt_adjusted_cox_est) %>% 
    #   summary(times = 5)
    
    
    ### STEP 4: Censoring weights ###
    
    expanded_df <- expandRows(df, "trt_weight")
    expanded_df <- expanded_df %>% arrange(obs_time, cens_ind)
    
    df2 <- expanded_df 
    
    expanded_df_trt1 <- expanded_df %>% filter(trt == 1)
    ipcw_expanded_df_trt1 <- ipcw(
      Hist(obs_time, 1-cens_ind) ~ x1 + x2 + x3,
      data = expanded_df_trt1,
      method = "cox",
      times = 5,
      subject.times = expanded_df_trt1$obs_time,
      keep = c("fit", "IPCW.times")
      )
    
    expanded_df_trt0 <- expanded_df %>% filter(trt == 0)
    ipcw_expanded_df_trt0 <- ipcw( 
      Hist(obs_time, 1-cens_ind) ~ x1 + x2 + x3,
      data = expanded_df_trt0,
      method = "cox",
      times = 5,
      subject.times = expanded_df_trt0$obs_time,
      keep = c("fit", "IPCW.times")
    )
    
    fitted_ipcw_trt1 <- ipcw_expanded_df_trt1$IPCW.times
    expanded_df_trt1$cens_weights = fitted_ipcw_trt1
    fitted_ipcw_trt0 <- ipcw_expanded_df_trt0$IPCW.times
    expanded_df_trt0$cens_weights = fitted_ipcw_trt0
    
    expanded_df <- bind_rows(expanded_df_trt1, expanded_df_trt0)
    
    weighted_cens_df <- 
      expanded_df %>%
      filter(cens_weights > 0)
    
    df3 <- weighted_cens_df
      
    ### STEP 5: Censoring-weighted Estimates ###
    trt_cens_adjusted_km_est <- 
      survfit(
        Surv(obs_time, 1-cens_ind) ~ trt, 
        data = weighted_cens_df,
        weights = cens_weights)
    
    cur_km_est_trt_cens_adj <- 
      trt_cens_adjusted_km_est %>% 
      summary(times = 5)
    # 
    # trt_cens_adjusted_cox_est <-
    #   coxph(
    #     Surv(obs_time, 1 - cens_ind) ~ strata(trt) + x1 + x2 + x3,
    #     data = weighted_cens_df,
    #     weights = cens_weights
    #   )
    # 
    # cur_cox_est_trt_cens_adj <- 
    #   survfit(trt_cens_adjusted_cox_est) %>% 
    #   summary(times = 5)
    # 
    
    #####################################
    #                                   #
    # Censoring  first, then treatment  #
    #                                   #
    #####################################
   
    ### STEP 1: IPCW ###
    df_trt1 <- df %>% filter(trt == 1)
    ipcw_df_trt1 <- ipcw(
      Hist(obs_time, 1-cens_ind) ~ x1 + x2 + x3,
      data = df_trt1,
      method = "cox",
      times = 5,
      subject.times = df_trt1$obs_time,
      keep = c("fit", "IPCW.times")
    )
    
    df_trt0 <- df %>% filter(trt == 0)
    ipcw_df_trt0 <- ipcw(
      Hist(obs_time, 1-cens_ind) ~ x1 + x2 + x3,
      data = df_trt0,
      method = "cox",
      times = 5,
      subject.times = df_trt0$obs_time,
      keep = c("fit", "IPCW.times")
    )
    
    fitted_ipcw_trt1 <- ipcw_df_trt1$IPCW.times
    df_trt1$cens_weights = fitted_ipcw_trt1
    fitted_ipcw_trt0 <- ipcw_df_trt0$IPCW.times
    df_trt0$cens_weights = fitted_ipcw_trt0
    
    df <- bind_rows(df_trt1, df_trt0)
    
    cens_df <- 
      df %>%
      filter(cens_weights > 0)
    
    df4 <- cens_df
    
    ### STEP 2: IPCW ESTIMATES ####
    
    cens_adjusted_km_est <- 
      survfit(
        Surv(obs_time, 1-cens_ind) ~ trt, 
        data =  cens_df,
        weights = cens_weights)
    
    
    cur_km_est_cens_adj <- 
      cens_adjusted_km_est %>% 
      summary(times = 5)
    
    # 
    # cens_adjusted_cox_est <-
    #   coxph(
    #     Surv(obs_time, 1 - cens_ind) ~ strata(trt) + x1 + x2 + x3,
    #     data = cens_df,
    #     weights = cens_weights
    #   )
    # 
    # cur_cox_est_cens_adj <- 
    #   survfit(cens_adjusted_cox_est) %>% 
    #   summary(times = 5)
    # 
    ### STEP 3: IPTW #########
    
    iptw_reg <- 
      glm(
        trt ~ x1 + x2 + x3, 
        family = "quasibinomial", 
        data = cens_df,
        weights = cens_weights)
    
    cens_df$prop_score <- iptw_reg$fitted.values
    
    cens_df <- 
      cens_df %>% mutate(
        trt_weight = 
          ifelse(trt == 1, 
                 1/prop_score, 
                 1/(1 - prop_score))
      )
    
    df5 <- cens_df
    
    ### STEP 5: Censoring-weighted Estimates ###
    
    cens_trt_adjusted_km_est <- 
      survfit(
        Surv(obs_time, 1-cens_ind) ~ trt, 
        data = cens_df,
        weights = trt_weight)
    
    cur_km_est_cens_trt_adj <- 
      cens_trt_adjusted_km_est %>% 
      summary(times = 5)
    
    # cens_trt_adjusted_cox_est <-
    #   coxph(
    #     Surv(obs_time, 1 - cens_ind) ~ strata(trt) + x1 + x2 + x3,
    #     data = cens_df,
    #     weights = trt_weight
    #   )
    # 
    # cur_cox_est_cens_trt_adj <- 
    #   survfit(cens_trt_adjusted_cox_est) %>% 
    #   summary(times = 5)
    # 
    
    
    ### trt adjustment all together ###
    
    #unadj_km_trteff_est <- cur_km_est_unadj$surv[2] - cur_km_est_unadj$surv[1]
    #unadj_cox_trteff_est <- cur_cox_est_unadj$surv[2] - cur_cox_est_unadj$surv[1]
    
    adj_trt_km_trteff_est <- cur_km_est_trt_adj$surv[2] - cur_km_est_trt_adj$surv[1]
    #adj_trt_cox_trteff_est <- cur_cox_est_trt_adj$surv[2] - cur_cox_est_trt_adj$surv[1]
    
    adj_cens_km_trteff_est <- cur_km_est_cens_adj$surv[2] - cur_km_est_cens_adj$surv[1]
    #adj_cens_cox_trteff_est <- cur_cox_est_cens_adj$surv[2] - cur_cox_est_cens_adj$surv[1]
    
    adj_trt_cens_km_trteff_est <- cur_km_est_trt_cens_adj$surv[2] - cur_km_est_trt_cens_adj$surv[1]
    #adj_trt_cens_cox_trteff_est <- cur_cox_est_trt_cens_adj$surv[2] - cur_cox_est_trt_cens_adj$surv[1]
    
    adj_cens_trt_km_trteff_est <- cur_km_est_cens_trt_adj$surv[2] - cur_km_est_cens_trt_adj$surv[1]
    #adj_cens_trt_cox_trteff_est <- cur_cox_est_cens_trt_adj$surv[2] - cur_cox_est_cens_trt_adj$surv[1]
    
    res_tib <-
      tibble(
        #est_order = 1:10,
        est_order = 1:4,
      est_name = 
        c("trt-adj KM", "cens-adj KM", "trt-cens-adj KM", "cens-trt-adj KM"),
          #"unadj KM", "trt-adj KM", "cens-adj KM", "trt-cens-adj KM", "cens-trt-adj KM"),
          #"unadj Cox", "trt-adj Cox", "cens-adj Cox", "trt-cens-adj Cox", "cens-trt-adj Cox"),
      est = c(
        #unadj_km_trteff_est,
              adj_trt_km_trteff_est,
              adj_cens_km_trteff_est,
              adj_trt_cens_km_trteff_est,
              adj_cens_trt_km_trteff_est
              
              # unadj_cox_trteff_est,
              # adj_trt_cox_trteff_est,
              # adj_cens_cox_trteff_est,
              # adj_trt_cens_cox_trteff_est,
              # adj_cens_trt_cox_trteff_est
              )
      )
      dat_tib <-
        tibble(
          df1 = df1 %>% list(),
          df2 = df2 %>% list(),
          df3 = df3 %>% list(),
          df4 = df4 %>% list(),
          df5 = df5 %>% list()
      )
    return(
      list(
        data = dat_tib,
        results = res_tib
      )
      )
 }

 samp_function <-
   function(data, id){
     samp_tib <- 
       sample_n(data, nrow(data), replace=TRUE) %>%
       mutate(
         boot_id = id
       )
     return(samp_tib)
   }
 
 extract_list_element <- function(list, i){
     return(list[[i]])
 }
 
bootstrapped_sim <- function(df, boot_n, true_param){
  # cur_tib <- gen_data(1000, 1, high_b, high_a, low_high_g)
  # df <- cur_tib[[1]]
  # boot_n = 5
  # true_param = 0.127
  df_tib <- 
    tibble(
      booted_index = 1:boot_n,
      df_list = map(booted_index, ~samp_function(id = ., data = df))
      )
  
  booted_tib <-
    df_tib %>%
    mutate(
      all_res = map(df_list, gen_res),
      data = map(all_res, ~extract_list_element(., 1)),
      res = map(all_res, ~extract_list_element(., 2))
    ) %>%
    select(data, res)
  
  out_tib <-  
    booted_tib  %>%
    mutate(
      res_orig = res
    ) %>%
    unnest(res) %>%
    group_by(est_name) %>%
    reframe(
      mean_est = mean(est, na.rm=TRUE), 
      se = sd(est, na.rm=TRUE),
      lower_ci = quantile(est, 0.025, na.rm=TRUE),
      upper_ci = quantile(est, 0.975, na.rm=TRUE),
      cover_ind = ifelse(true_param <= upper_ci & true_param >= lower_ci, 1, 0),
      data = data,
      res = res_orig
    )
  
  return(out_tib)
}

n_sims <- 100
n_rows <- 1000

#scenario 1 
run_sim <- function(
    te, beta_vec, alpha_vec, gamma_vec, n_s = n_sims, n_r = n_rows
  ){
  # te <- 0.12
  # beta_vec = high_b
  # alpha_vec = high_a
  # gamma_vec = low_high_g

    tibs <- 
      gen_data(n_rows, 
               n_sims, 
               betas = beta_vec, 
               alphas = alpha_vec, 
               gammas = gamma_vec)
      
    sim <- 
        tibble(
        simmed_tibs = tibs,
        sim_res = map(simmed_tibs, ~bootstrapped_sim(df = ., boot_n = n_sims, true_param = te), .progress=TRUE)
      ) %>%
      select(sim_res) %>%
      unnest(sim_res)
    
    sim_dat <-
      sim %>% 
      select(data, res) %>%
      unique()
    
    sim_res_sum <- 
      sim %>% 
      select(est_name, mean_est, se, lower_ci, upper_ci, cover_ind) %>% 
      unique() %>%
      group_by(est_name) %>%
      summarise(
        avg_te = mean(mean_est, na.rm=TRUE),
        avg_se = mean(se, na.rm=TRUE),
        coverage_prop = mean(cover_ind, na.rm=TRUE),
        relative_bias = mean(mean_est - te, na.rm=TRUE)/te
      )
    beep()
    return(list(sim_dat, sim_res_sum))
}

scen1_sum <- run_sim(0.12, high_b, high_a, low_high_g)
scen1_sum[[2]] <- scen1_sum[[2]] %>% mutate(scenario = "scen1")
save(scen1_sum, file="scen1_sum_100.RData")
rm(scen1_sum)

scen2_sum <- run_sim(0.12, med_b, high_a, low_high_g) 
scen2_sum[[2]] <- scen2_sum[[2]] %>% mutate(scenario = "scen2")
save(scen2_sum, file="scen2_sum_100.RData")
rm(scen2_sum)

scen3_sum <- run_sim(0.12, low_b, high_a, low_high_g) 
scen3_sum[[2]] <- scen3_sum[[2]] %>% mutate(scenario = "scen3")
save(scen3_sum, file="scen3_sum_100.RData")
rm(scen3_sum)

scen4_sum <- run_sim(0.014, high_b, low_a, low_low_g) 
scen4_sum[[2]] <- scen4_sum[[2]] %>% mutate(scenario = "scen4")
save(scen4_sum, file="scen4_sum_100.RData")
rm(scen4_sum)

scen5_sum <- run_sim(0.014, med_b, low_a, low_low_g)
scen5_sum[[2]] <- scen5_sum[[2]] %>% mutate(scenario = "scen5")
save(scen5_sum, file="scen5_sum_100.RData")
rm(scen5_sum)

scen6_sum <- run_sim(0.014, low_b, low_a, low_low_g)
scen6_sum[[2]] <- scen6_sum[[2]] %>% mutate(scenario = "scen6")
save(scen6_sum, file="scen6_sum_100.RData")
rm(scen6_sum)

scen7_sum <- run_sim(0.12, high_b, high_a, high_high_g) 
scen7_sum[[2]] <- scen7_sum[[2]] %>% mutate(scenario = "scen7")
save(scen7_sum, file="scen7_sum_100.RData")
rm(scen7_sum)

scen8_sum <- run_sim(0.12, med_b, high_a, high_high_g) 
scen8_sum[[2]] <- scen8_sum[[2]] %>% mutate(scenario = "scen8")
save(scen8_sum, file="scen8_sum_100.RData")
rm(scen8_sum)

scen9_sum <- run_sim(0.12, low_b, high_a, high_high_g) 
scen9_sum[[2]] <- scen9_sum[[2]] %>% mutate(scenario = "scen9")
save(scen9_sum, file="scen9_sum_100.RData")
rm(scen9_sum)

scen10_sum <- run_sim(0.014, high_b, low_a, high_low_g) 
scen10_sum[[2]] <- scen10_sum[[2]] %>% mutate(scenario = "scen10")
save(scen10_sum, file="scen10_sum_100.RData")
rm(scen10_sum)

scen11_sum <- run_sim(0.014, med_b, low_a, high_low_g)
scen11_sum[[2]] <- scen11_sum[[2]] %>% mutate(scenario = "scen11")
save(scen11_sum, file="scen11_sum_100.RData")
rm(scen11_sum)

scen12_sum <- run_sim(0.014, low_b, low_a, high_low_g)
scen12_sum[[2]] <- scen12_sum[[2]] %>% mutate(scenario = "scen12")
save(scen12_sum, file="scen12_sum_100.RData")
rm(scen12_sum)

beep(8)

# 
# all_scen_sum_100 <- bind_rows(
#   scen1_sum,
#   scen2_sum,
#   scen3_sum,
#   scen4_sum,
#   scen5_sum,
#   scen6_sum,
#   scen7_sum,
#   scen8_sum,
#   scen9_sum,
#   scen10_sum,
#   scen11_sum,
#   scen12_sum) %>%
#   mutate(
#     avg_te = round(avg_te, digits = 3)
#   )
# #save(all_scen_sum_50, file="all_scen_sum_50.RData")
# 
# cur_sum_tib <- all_scen_sum_50 %>%
#   filter(est_name == "cens-trt-adj KM") %>%
#   select(scenario, avg_te, relative_bias, coverage_prop)
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = coverage_prop,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("KM",est_name)
#   ) %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = relative_bias,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("KM",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = avg_te,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("KM",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = avg_se,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("KM",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = coverage_prop,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("KM",est_name)
#   ) %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = relative_bias,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("KM",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = avg_te,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("KM",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = avg_se,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("KM",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)

# 
# 
# all_scen_sum_10 <- bind_rows(
#   scen1_sum,
#   scen2_sum,
#   scen3_sum,
#   scen4_sum,
#   scen5_sum,
#   scen6_sum,
#   scen7_sum,
#   scen8_sum,
#   scen9_sum,
#   scen10_sum,
#   scen11_sum,
#   scen12_sum) 
# 
# 
# 
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = coverage_prop,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("Cox",est_name)
#   ) %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = relative_bias,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("Cox",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = avg_te,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("Cox",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = avg_se,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("Cox",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = coverage_prop,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("Cox",est_name)
#   ) %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = relative_bias,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("Cox",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = avg_te,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("Cox",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
# 
# all_scen_sum_50 %>%
#   pivot_wider(
#     values_from = avg_se,
#     names_from = scenario,
#     id_cols = est_name
#   ) %>%
#   filter(
#     grepl("Cox",est_name)
#   )  %>%
#   bind_cols(
#     order = c(4,5,2,3,1)
#   ) %>%
#   arrange(order) %>%
#   dplyr::select(-order)
