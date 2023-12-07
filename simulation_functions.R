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

# example run
#
 cur_tib <- gen_data(1000, 
                   10, 
                    betas = high_b, 
                    alphas = high_a, 
                    gammas = high_high_g)

ipw <- function(df){
  
 df <- cur_tib[[2]]
  df %>% 
    summarise(
      censoring_rate = sum(cens_ind)/n()
    )
  
  trt_res <- glm(trt ~ x1 + x2 + x3, family = "binomial", data = df)
  cens_res <- glm(cens_ind ~ x1 + x2 + x3, family = "binomial", data = df)
  
  summary(cens_res)
  summary(trt_res)
  df$trt_weights <- 
    ifelse(df$trt == 1, 
           1/trt_res$fitted.values, 
           1/(1 - trt_res$fitted.values))
  df$cens_weights <- 
    ifelse(df$cens_ind == 1, 
           1/cens_res$fitted.values, 
           1/(1 - cens_res$fitted.values))
  
  df %>%
    group_by(trt) %>% 
    summarise(
      prop_x1 = mean(x1),
      weighted_prop_x1 = mean(x1 * trt_weights),
      mean_x2 = mean(x2),
      weighted_mean_x2 = mean(x2*trt_weights),
      mean_x3 = mean(x3),
      weighted_mean_x3 = mean(x3*trt_weights)
  )
  
  df %>%
    filter(cens_weights < quantile(cens_weights, 0.99) & cens_weights > quantile(cens_weights, 0.01)) %>%
    group_by(cens_ind) %>% 
    summarise(
      prop_x1 = mean(x1),
      weighted_prop_x1 = mean(x1*cens_weights),
      mean_x2 = mean(x2),
      weighted_mean_x2 = mean(x2*cens_weights),
      mean_x3 = mean(x3),
      weighted_mean_x3 = mean(x3*cens_weights)
    )
}

tibble(
  simmed_tibs = cur_tib,
  res = map(simmed_tibs, nrow) %>% unlist()
)









