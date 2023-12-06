# ####################################
# 
# Program: propensity_score_functions.R
#
# Description: Contains functions for propensity score generation
# 
# ####################################

# ###################################################
# gen_prop_score_beta0()
# - used to find beta0 given overlap parameters and 
#   desired treatment proportion
# - pass in covariate list and model parameters
# - returns list containing selected beta0, full search
#   procedure results, and ggplot of beta0
# ###################################################

gen_prop_score_beta0 <- function(n_prop,
                                 p=0.55, 
                                 b, 
                                 thresh = 0.01){

  beta1 <- b[[1]]
  beta2 <- b[[2]]
  
  lower_bound <- -10
  upper_bound <- 10
  step <- 0.01
  
  cur_b0_step <- lower_bound
  
  cur_res <- 
    tibble(
    )
  
  pb <- progress_bar$new(
    total = (upper_bound - lower_bound)/step + 1, 
    format = "running sim of size :total... [:bar]   :percent completed; eta: :eta"
  )
  while(cur_b0_step <= upper_bound){
    pb$tick()
    x1 <- gen_cur_covs(n_prop)
    
    x1 <- x[[1]]
    x2 <- x[[2]]
    x3 <- x[[3]]
    
    cur_b0 <- cur_b0_step
    cur_b1 <- beta1
    cur_b2 <- beta2
    
    cur_to_expit <- cur_b0 + cur_b1 * x1 + cur_b2 * x2
    
    probs <- exp(cur_to_expit)/(1 + exp(cur_to_expit))
    cur_probs <- sum(probs > 0.5)/sum(!is.na(probs))
    
    cur_res <- bind_rows(cur_res, tibble(
      b0 = cur_b0,
      b1 = cur_b1,
      b2 = cur_b2,
      probs = cur_probs
    ))
    
    cur_b0_step <- cur_b0_step + step
  }
  b0_fin <- 
    cur_res %>% 
    filter(probs <= p + thresh & probs >= p - thresh) %>% 
    pull(b0) %>% 
    median()
  
  cur_res_gg <-
    cur_res %>% 
    ggplot(
      aes(x = b0, y = probs)
    ) + geom_line() +
    geom_hline(yintercept = p, col = "red") +
    labs(
      x = "Beta_0",
      y = "Probability of Treatment",
      title = "Estimated Treatment Probability by Beta_0"
    )
  
  return(list(b0_emp = b0_fin, cur_res = cur_res, cur_res_gg = cur_res_gg))
}

# ### parameter setting beta 0 ###

#high_b <- list(0.1, 0.1)
#med_b <- list(1, 1)
#low_b <- list(3, 3)
#cur_n <- 100000

# ### low ###
#beta_0_low <- gen_prop_score_beta0(cur_n, b=low_b)

# ### med ###
#beta_0_med <- gen_prop_score_beta0(cur_n, b=med_b)

# ### high ###
#beta_0_high <- gen_prop_score_beta0(cur_n, b=high_b)

# ###################################################
# gen_prop_scores()
# - used to append propensity scores to covariate
#   dataframe
# - pass in beta list and covariate dataframe
# - returns covariate dataframe with prop scores and
#   treatment assignment appended
# ###################################################

gen_prop_scores <- function(b, covs = cov_df){
  
  b0 <- b[[1]]
  b1 <- b[[2]]
  b2 <- b[[3]]
  n <- nrow(covs)
  
  out_df <- 
    covs %>%
    mutate(
      prop_score = 
        exp(b0 + b1*x1 + b2*x2)/
        (1 + exp(b0 + b1*x1 + b2*x2))
    )
  
  out_df$trt <- rbinom(n, 1, out_df$prop_score)
  
  return(out_df)
}
