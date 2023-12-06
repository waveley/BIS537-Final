# ####################################
# 
# Program: propensity_score_fns.R
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

gen_prop_score_beta0 <- function(x=xlist, 
                                 p=0.55, 
                                 b, 
                                 thresh = 0.01){
  x1 <- x[[1]]
  x2 <- x[[2]]
  x3 <- x[[3]]
  beta1 <- b[[1]]
  beta2 <- b[[2]]
  beta3 <- b[[3]]
  
  lower_bound <- -10
  upper_bound <- 10
  step <- 0.01
  
  cur_b0_step <- lower_bound
  
  cur_res <- 
    tibble(
    )
  
  while(cur_b0_step <= upper_bound){
    
    x1 <- 
    
    cur_b0 <- cur_b0_step
    cur_b1 <- beta1
    cur_b2 <- beta2
    cur_b3 <- beta3
    
    cur_to_expit <- cur_b0 + cur_b1 * x1 + cur_b2 * x2
    
    probs <- exp(cur_to_expit)/(1 + exp(cur_to_expit))
    cur_probs <- sum(probs > 0.5)/sum(!is.na(probs))
    
    cur_res <- bind_rows(cur_res, tibble(
      b0 = cur_b0,
      b1 = cur_b1,
      b2 = cur_b2,
      b3 = cur_b3,
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
