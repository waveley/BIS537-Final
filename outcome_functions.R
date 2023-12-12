# ################################################
# 
# Program: outcome_functions.R
#
# Description: Contains functions for outcome 
#              and censoring generation
# 
# ################################################


# ###################################################
# gen_outcomes()
# - pass in covariate dataframe and model parameters
# - returns covariate dataframe with outcomes
# ###################################################

gen_outcomes <- 
  function(covs, 
           lambda=0.0001, 
           nu=3, 
           a = list(a0, a1, a3)){
  a0 <- a[[1]]
  a1 <- a[[2]]
  a3 <- a[[3]]
  n <- nrow(covs)
  
  out_df <- 
    covs %>% 
    bind_cols(tibble(u_surv = runif(n))) %>%
    mutate(L = a0*trt + a1*x1 + a3*x3,
           surv_time = (-log(u_surv)/(lambda * exp(L)))^(1/nu)
    ) %>% 
    dplyr::select(
      -L, -u_surv
    )
  return(out_df)
}

# ###################################################
# gen_censor()
# - pass in covariate dataframe and model parameters
#   and survival time
# - returns covariate dataframe with censoring time
#   and censoring indicator
# ###################################################

gen_censor <- function(covs, 
                       lambda=0.0001, 
                       nu=3, 
                       g = list(g0, g2, g3)){
  g0 <- g[[1]]
  g2 <- g[[2]]
  g3 <- g[[3]]
  n <- nrow(covs)
  
  out_df <- 
    covs %>% 
    bind_cols(tibble(u_cens = runif(n))) %>%
    mutate(K = g0 + g2*x2 + g3*x3,
           cens_time = (-log(u_cens)/(lambda*exp(K))),
           cens_ind = ifelse(cens_time < surv_time, 1, 0)
    ) %>% 
    dplyr::select(
      -K, -u_cens
    )
  return(out_df)
}

# ###################################################
# gen_censor_gamma0()
# - used to find gamma0 given parameters, desired censoring
#   proportion, covariates, and true survival time
#   desired treatment proportion
# - pass in covariate list and model parameters
# - returns list containing selected beta0, full search
#   procedure results, and ggplot of beta0
# ###################################################

gen_censor_gamma0 <- function(n_gen=10000, 
                              p_cens=0.25, 
                              g, 
                              b,
                              a,
                              lambda = 0.0001, 
                              nu = 3,
                              thresh = 0.01){
  lower_bound <- -10
  upper_bound <- 10
  step <- 0.01
  
  cur_g0_step <- lower_bound
  
  cur_res <- 
    tibble(
    )
  
  pb <- progress_bar$new(
    total = (upper_bound - lower_bound)/step + 1, 
    format = "running sim of size :total... [:bar]   :percent completed; eta: :eta"
    )
  
  while(cur_g0_step <= upper_bound){
    pb$tick()
    x <- gen_cur_covs(n_gen)
    x1 <- x[[1]]
    x2 <- x[[2]]
    x3 <- x[[3]]
    
    gamma2 <- g[[1]]
    gamma3 <- g[[2]]
    
    cur_g0 <- cur_g0_step
    cur_g2 <- gamma2
    cur_g3 <- gamma3
    cur_g <- list(cur_g0, cur_g2, cur_g3)
      
    # generate censoring time
    cur_cov <- gen_prop_scores(b, tibble(x1, x2, x3))
    cur_out_df <- gen_outcomes(cur_cov, a = a)
    cur_out_cens_df <- gen_censor(cur_out_df, g = cur_g)
    
    cur_cens_prob <- sum(cur_out_cens_df %>% pull(cens_ind))/n_gen
    
    #print(cur_cens_prob)
    cur_res <- bind_rows(cur_res, tibble(g0 = cur_g0_step, prop_cens = cur_cens_prob))
    cur_g0_step <- cur_g0_step + step
  }
  selected_g0 <- cur_res %>%
    filter(prop_cens <= p_cens + thresh & prop_cens >= p_cens - thresh) %>%
    arrange(-g0) %>%
    slice(1)
  
  return(list(cur_res, selected_g0))
}

### parameter setting -- gamma0 ###
# 
# g <- list(2, 4)

# # scenario 1
# 
# cur_gamma_sim_scen1 <- 
#   gen_censor_gamma0(
#     g = list(2,4),
#     b = high_b,
#     a = high_a,
#     p_cens = low_p_cens
#     )
# 
# # scenario 2
# 
# cur_gamma_sim_scen2 <- 
#   gen_censor_gamma0(
#     g = g,
#     b = med_b,
#     a = high_a,
#     p_cens = low_p_cens
#   )
# 
# # scenario 3
# 
# cur_gamma_sim_scen3 <- 
#   gen_censor_gamma0(
#     g = g,
#     b = low_b,
#     a = high_a,
#     p_cens = low_p_cens
#   )
# 
# # scenario 4
# 
# cur_gamma_sim_scen4 <- 
#   gen_censor_gamma0(
#     g = g,
#     b = high_b,
#     a = low_a,
#     p_cens = low_p_cens
#   )
# 
# # scenario 5
# 
# cur_gamma_sim_scen5 <- 
#   gen_censor_gamma0(
#     g = g,
#     b = med_b,
#     a = low_a,
#     p_cens = low_p_cens
#   )
# 
# # scenario 6
# 
# cur_gamma_sim_scen6 <- 
#   gen_censor_gamma0(
#     g = g,
#     b = low_b,
#     a = low_a,
#     p_cens = low_p_cens
#   )
# 
# # scenario 7
# 
# cur_gamma_sim_scen7 <- 
#   gen_censor_gamma0(
#     g = g,
#     b = high_b,
#     a = high_a,
#     p_cens = high_p_cens
#   )
# 
# # scenario 8
# 
# cur_gamma_sim_scen8 <- 
#   gen_censor_gamma0(
#     g = g,
#     b = med_b,
#     a = high_a,
#     p_cens = high_p_cens
#   )
# 
# # scenario 9
# 
# cur_gamma_sim_scen9 <- 
#   gen_censor_gamma0(
#     g = g,
#     b = low_b,
#     a = high_a,
#     p_cens = high_p_cens
#   )
# 
# # scenario 10
# 
# cur_gamma_sim_scen10 <- 
#   gen_censor_gamma0(
#     g = g,
#     b = high_b,
#     a = low_a,
#     p_cens = high_p_cens
#   )
# 
# # scenario 11
# 
# cur_gamma_sim_scen11 <- 
#   gen_censor_gamma0(
#     g = g,
#     b = med_b,
#     a = low_a,
#     p_cens = high_p_cens
#   )
# 
# # scenario 12
# 
# cur_gamma_sim_scen12 <- 
#   gen_censor_gamma0(
#     g = g,
#     b = low_b,
#     a = low_a,
#     p_cens = high_p_cens
#   )


# ### true treatment effects ###
# 
# library(progress)
# library(beepr)
# library(patchwork)
# 
# run_true_survival <- function(n, params, t_end = 5){
#   
#   alpha_0 <- params[[1]]
#   alpha_1 <- params[[2]]
#   alpha_3 <- params[[3]]
#   
#   lambda <- params[[4]]
#   v <- params[[5]]
  
  ###### define confounders here #####
  # x <- gen_cur_covs(n)
  # x_1 <- x[[1]]
  # x_3 <- x[[3]]
  
#   ###### define l here #####
#   l <- alpha_0 + alpha_1*x_1 + alpha_3*x_3 
#   
#   t <- c()
#   
#   pb <- progress_bar$new(total = n, format = "running sim of size :total... [:bar]   :percent completed; eta: :eta")
#   
#     for(i in 1:n){
#       pb$tick()
#       cur_samp <- runif(n)
#       calculated_samp <- (-log(cur_samp)/(lambda*exp(l)))^(1/v)
#       t <- mean(calculated_samp > t_end)
#     }
#   
#     t_hat <- mean(t)
#     beep()
#   return(t_hat)
# }
# 
# m <- 100000
# 
# gen_expected_ace <- function(n, params1, params0){
#   t_1_true_val <- run_true_survival(n, params = params1) %>% round(digits = 3)
#   t_0_true_val <- run_true_survival(n, params = params0) %>% round(digits = 3)
#   return(t_1_true_val - t_0_true_val)
# }
# 
# high_te <- gen_expected_ace(n = m, list(-1, 2, 3, 0.0001, 3), list(0, 2, 3, 0.0001, 3))
# low_te <- gen_expected_ace(n = m, list(-0.1, 2, 3, 0.0001, 3), list(0, 2, 3, 0.0001, 3))
# high_te
# low_te
