# ####################################
# 
# Program: outcome_functions.R
#
# Description: Contains functions for outcome 
#              and censoring generation
# 
# ####################################


# ###################################################
# gen_outcomes()
# - pass in covariate dataframe and model parameters
# - returns covariate dataframe with outcomes
# ###################################################

gen_outcomes <- function(covs, lambda=0.0001, nu=3, a = list(a0, a1, a3)){
  a0 <- a[[1]]
  a1 <- a[[2]]
  a3 <- a[[3]]
  n <- nrow(covs)
  
  out_df <- 
    covs %>% 
    bind_cols(tibble(u_surv = runif(n))) %>%
    mutate(L = a0*trt + a1*x1 + a3*x3,
           surv_time = (-log(u_surv)/(lambda * exp(L)))^(1/nu)
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

gen_censor <- function(covs, lambda=0.0001, nu=3, g = list(g0, g2, g3)){
  g0 <- g[[1]]
  g2 <- g[[2]]
  g3 <- g[[3]]
  n <- nrow(covs)
  
  out_df <- 
    covs %>% 
    bind_cols(tibble(u_cens = runif(n))) %>%
    mutate(K = exp(g0 + g2*x2 + g3*x3),
           cens_time = (-log(u_cens)/(lambda*exp(K))),
           cens_ind = ifelse(cens_time < surv_time, 1, 0),
           obs_time = min(surv_time, cens_time)
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

gen_censor_gamma0 <- function(n_gen=100000, 
                              p_cens=0.25, 
                              g, 
                              b,
                              a,
                              lambda = 0.0001, 
                              nu = 3,
                              thresh = 0.01){
  lower_bound <- -10
  upper_bound <- 1
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

# scenario 1
cur_gamma_sim_scen1 <- gen_censor_gamma0(g = list(2, 4),
                                   b = list(0.141, 0.1, 0.1),
                                   a = list(-1, 2, 3))

# scenario 2
cur_gamma_sim_scen2- gen_censor_gamma0(g = list(2, 4),
                                   b = list(-0.399, 1, 1),
                                   a = list(-1, 2, 3))

# scenario 3
cur_gamma_sim_scen3 <- gen_censor_gamma0(g = list(2, 4),
                                       b = list(-1.598, 1, 1),
                                       a = list(-1, 2, 3))
# scenario 4
cur_gamma_sim_scen4 <- gen_censor_gamma0(g = list(2, 4),
                                         b = list(0.141, 0.1, 0.1),
                                         a = list(-1, 2, 3),
                                         p_cens = 0.5)

temp_df <- cur_gamma_sim[[1]]
