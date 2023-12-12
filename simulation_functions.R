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
   unadjusted_km_est <- 
      survfit(
        Surv(obs_time, 1-cens_ind) ~ trt, 
        data = df)
    
    cur_km_est_unadj <- 
      unadjusted_km_est %>% 
      summary(times = 5)
    
    unadjusted_cox_est <-
      coxph(
        Surv(obs_time, 1 - cens_ind) ~ strata(trt) + x1 + x2 + x3,
        data = df
      )
    
    cur_cox_est_unadj <- 
      survfit(unadjusted_cox_est) %>% 
      summary(times = 5)
    
    cur_cox_est_unadj
    cur_km_est_unadj
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
  
    ### STEP 3: Treatment-weighted Estimates ###
    
    trt_adjusted_km_est <- 
      survfit(
        Surv(obs_time, 1-cens_ind) ~ trt, 
        data = df,
        weights = trt_weight)
    
    
    cur_km_est_trt_adj <- 
      trt_adjusted_km_est %>% 
      summary(times = 5)
     
    
    trt_adjusted_cox_est <-
      coxph(
        Surv(obs_time, 1 - cens_ind) ~ strata(trt) + x1 + x2 + x3,
        data = df,
        weights = trt_weight
      )
    
    cur_cox_est_trt_adj <- 
      survfit(trt_adjusted_cox_est) %>% 
      summary(times = 5)
    
    
    ### STEP 4: Censoring weights ###
    
    expanded_df <- expandRows(df, "trt_weight")
    expanded_df <- expanded_df %>% arrange(obs_time, cens_ind)
    
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
      
    ### STEP 5: Censoring-weighted Estimates ###
    trt_cens_adjusted_km_est <- 
      survfit(
        Surv(obs_time, 1-cens_ind) ~ trt, 
        data = weighted_cens_df,
        weights = cens_weights)
    
    cur_km_est_trt_cens_adj <- 
      trt_cens_adjusted_km_est %>% 
      summary(times = 5)
    
    trt_cens_adjusted_cox_est <-
      coxph(
        Surv(obs_time, 1 - cens_ind) ~ strata(trt) + x1 + x2 + x3,
        data = weighted_cens_df,
        weights = cens_weights
      )
    
    cur_cox_est_trt_cens_adj <- 
      survfit(trt_cens_adjusted_cox_est) %>% 
      summary(times = 5)
    
    
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
    
    ### STEP 2: IPCW ESTIMATES ####
    
    cens_adjusted_km_est <- 
      survfit(
        Surv(obs_time, 1-cens_ind) ~ trt, 
        data =  cens_df,
        weights = cens_weights)
    
    
    cur_km_est_cens_adj <- 
      cens_adjusted_km_est %>% 
      summary(times = 5)
    
    
    cens_adjusted_cox_est <-
      coxph(
        Surv(obs_time, 1 - cens_ind) ~ strata(trt) + x1 + x2 + x3,
        data = cens_df,
        weights = cens_weights
      )
    
    cur_cox_est_cens_adj <- 
      survfit(cens_adjusted_cox_est) %>% 
      summary(times = 5)
    
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
    
    ### STEP 5: Censoring-weighted Estimates ###
    
    cens_trt_adjusted_km_est <- 
      survfit(
        Surv(obs_time, 1-cens_ind) ~ trt, 
        data = cens_df,
        weights = trt_weight)
    
    cur_km_est_cens_trt_adj <- 
      cens_trt_adjusted_km_est %>% 
      summary(times = 5)
    
    cens_trt_adjusted_cox_est <-
      coxph(
        Surv(obs_time, 1 - cens_ind) ~ strata(trt) + x1 + x2 + x3,
        data = cens_df,
        weights = trt_weight
      )
    
    cur_cox_est_cens_trt_adj <- 
      survfit(cens_trt_adjusted_cox_est) %>% 
      summary(times = 5)
    
    
    
    ### trt adjustment all together ###
    
    unadj_km_trteff_est <- cur_km_est_unadj$surv[2] - cur_km_est_unadj$surv[1]
    unadj_cox_trteff_est <- cur_cox_est_unadj$surv[2] - cur_cox_est_unadj$surv[1]
    
    adj_trt_km_trteff_est <- cur_km_est_trt_adj$surv[2] - cur_km_est_trt_adj$surv[1]
    adj_trt_cox_trteff_est <- cur_cox_est_trt_adj$surv[2] - cur_cox_est_trt_adj$surv[1]
    
    adj_cens_km_trteff_est <- cur_km_est_cens_adj$surv[2] - cur_km_est_cens_adj$surv[1]
    adj_cens_cox_trteff_est <- cur_cox_est_cens_adj$surv[2] - cur_cox_est_cens_adj$surv[1]
    
    adj_trt_cens_km_trteff_est <- cur_km_est_trt_cens_adj$surv[2] - cur_km_est_trt_cens_adj$surv[1]
    adj_trt_cens_cox_trteff_est <- cur_cox_est_trt_cens_adj$surv[2] - cur_cox_est_trt_cens_adj$surv[1]
    
    adj_cens_trt_km_trteff_est <- cur_km_est_cens_trt_adj$surv[2] - cur_km_est_cens_trt_adj$surv[1]
    adj_cens_trt_cox_trteff_est <- cur_cox_est_cens_trt_adj$surv[2] - cur_cox_est_cens_trt_adj$surv[1]
    
    res_tib <-
      tibble(
        est_order = 1:10,
      est_name = 
        c("unadj KM", "trt-adj KM", "cens-adj KM", "trt-cens-adj KM", "cens-trt-adj KM",
          "unadj Cox", "trt-adj Cox", "cens-adj Cox", "trt-cens-adj Cox", "cens-trt-adj Cox"),
      est = c(unadj_km_trteff_est,
              adj_trt_km_trteff_est,
              adj_cens_km_trteff_est,
              adj_trt_cens_km_trteff_est,
              adj_cens_trt_km_trteff_est,
              
              unadj_cox_trteff_est,
              adj_trt_cox_trteff_est,
              adj_cens_cox_trteff_est,
              adj_trt_cens_cox_trteff_est,
              adj_cens_trt_cox_trteff_est
              )
      )
    return(res_tib)
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
 
bootstrapped_sim <- function(df, boot_n, true_param){
  
  df_tib <- 
    tibble(
      booted_index = 1:boot_n,
      df_list = map(booted_index, ~samp_function(id = ., data = df))
      )
  
  booted_tib <-
    df_tib %>%
    mutate(
      all_res = map(df_list, gen_res)
    ) %>%
    dplyr::select(-df_list)
  
  out_tib <-  
    booted_tib %>%
    unnest(all_res) %>% 
    group_by(est_name) %>%
    summarise(
      mean_est = mean(est), 
      se = sd(est),
      lower_ci = mean_est - 1.96*se,
      upper_ci = mean_est + 1.96*se,
      cover_ind = ifelse(true_param <= upper_ci & true_param >= lower_ci, 1, 0)
    )
  
  return(out_tib)
}

n_sims <- 100
n_rows <- 1000

#scenario 1 
run_sim <- function(
    te, beta_vec, alpha_vec, gamma_vec, n_s = n_sims, n_r = n_rows
  ){
    tibs <- 
      gen_data(n_rows, 
               n_sims, 
               betas = beta_vec, 
               alphas = alpha_vec, 
               gammas = gamma_vec)
      
    sim <- 
        tibble(
        simmed_tibs = tibs,
        res = map(simmed_tibs, ~bootstrapped_sim(df = ., boot_n = n_sims, true_param = te), .progress=TRUE)
      ) 
    
    cur_res_tib <- 
      sim %>% 
      unnest(res) %>% 
      dplyr::select(-simmed_tibs) %>%
      group_by(est_name) %>%
      summarise(
        avg_te = mean(mean_est),
        avg_se = mean(se),
        coverage_prop = mean(cover_ind),
        relative_bias = mean(mean_est - te)/te
      )
    
    return(cur_res_tib)
}

scen1_sum <- run_sim(0.12, high_b, high_a, low_high_g) %>% mutate(scenario = "scen1")
scen2_sum <- run_sim(0.12, med_b, high_a, low_high_g) %>% mutate(scenario = "scen2")
scen3_sum <- run_sim(0.12, low_b, high_a, low_high_g) %>% mutate(scenario = "scen3")
scen4_sum <- run_sim(0.014, high_b, low_a, low_low_g) %>% mutate(scenario = "scen4")
scen5_sum <- run_sim(0.014, med_b, low_a, low_low_g) %>% mutate(scenario = "scen5")
scen6_sum <- run_sim(0.014, low_b, low_a, low_low_g) %>% mutate(scenario = "scen6")
scen7_sum <- run_sim(0.12, high_b, high_a, high_high_g) %>% mutate(scenario = "scen7")
scen8_sum <- run_sim(0.12, med_b, high_a, high_high_g) %>% mutate(scenario = "scen8")
scen9_sum <- run_sim(0.12, low_b, high_a, high_high_g) %>% mutate(scenario = "scen9")
scen10_sum <- run_sim(0.014, high_b, low_a, high_low_g) %>% mutate(scenario = "scen10")
scen11_sum <- run_sim(0.014, med_b, low_a, high_low_g) %>% mutate(scenario = "scen11")
scen12_sum <- run_sim(0.014, low_b, low_a, high_low_g) %>% mutate(scenario = "scen12")

beep()

all_scen_sum_100<- bind_rows(
  scen1_sum,
  scen2_sum,
  scen3_sum,
  scen4_sum,
  scen5_sum,
  scen6_sum,
  scen7_sum,
  scen8_sum,
  scen9_sum,
  scen10_sum,
  scen11_sum,
  scen12_sum) %>%
  mutate(
    avg_te = round(avg_te, digits = 3)
  )
#save(all_scen_sum_50, file="all_scen_sum_50.RData")

cur_sum_tib <- all_scen_sum_50 %>%
  filter(est_name == "cens-trt-adj KM") %>%
  select(scenario, avg_te, relative_bias, coverage_prop)

all_scen_sum_50 %>%
  pivot_wider(
    values_from = coverage_prop,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("KM",est_name)
  ) %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)


all_scen_sum_50 %>%
  pivot_wider(
    values_from = relative_bias,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("KM",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)

all_scen_sum_50 %>%
  pivot_wider(
    values_from = avg_te,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("KM",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)

all_scen_sum_50 %>%
  pivot_wider(
    values_from = avg_se,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("KM",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)

all_scen_sum_50 %>%
  pivot_wider(
    values_from = coverage_prop,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("KM",est_name)
  ) %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)


all_scen_sum_50 %>%
  pivot_wider(
    values_from = relative_bias,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("KM",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)

all_scen_sum_50 %>%
  pivot_wider(
    values_from = avg_te,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("KM",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)

all_scen_sum_50 %>%
  pivot_wider(
    values_from = avg_se,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("KM",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)

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


all_scen_sum_50 %>%
  pivot_wider(
    values_from = coverage_prop,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("Cox",est_name)
  ) %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)


all_scen_sum_50 %>%
  pivot_wider(
    values_from = relative_bias,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("Cox",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)

all_scen_sum_50 %>%
  pivot_wider(
    values_from = avg_te,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("Cox",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)

all_scen_sum_50 %>%
  pivot_wider(
    values_from = avg_se,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("Cox",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)

all_scen_sum_50 %>%
  pivot_wider(
    values_from = coverage_prop,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("Cox",est_name)
  ) %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)


all_scen_sum_50 %>%
  pivot_wider(
    values_from = relative_bias,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("Cox",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)

all_scen_sum_50 %>%
  pivot_wider(
    values_from = avg_te,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("Cox",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)

all_scen_sum_50 %>%
  pivot_wider(
    values_from = avg_se,
    names_from = scenario,
    id_cols = est_name
  ) %>%
  filter(
    grepl("Cox",est_name)
  )  %>%
  bind_cols(
    order = c(4,5,2,3,1)
  ) %>%
  arrange(order) %>%
  dplyr::select(-order)
