
library(tidyverse)

gen_cur_covs <- function(n=1000){
  
  x1 <- runif(n, min = 0, max = 1)
  x1[x1 < 0.6] <- 1
  x1[x1 < 1] <- 0
  
  x2 <- rnorm(n, mean = 0, sd = 1)
  
  x3 <- rgamma(n, 1, 1)
  
  cov_df <- tibble(x1 = x1,
                   x2 = x2,
                   x3 = x3)
  
  xlist <-  list(x1, x2, x3)
  return(xlist)
}


gen_cur_covs_tib <- function(n=1000){
  
  x1 <- runif(n, min = 0, max = 1)
  x1[x1 < 0.6] <- 1
  x1[x1 < 1] <- 0
  
  x2 <- rnorm(n, mean = 0, sd = 1)
  
  x3 <- rgamma(n, 1, 1)
  
  cov_df <- tibble(x1 = x1,
                   x2 = x2,
                   x3 = x3)
  
  xlist <-  list(x1, x2, x3)
  return(cov_df)
}

# parameter setting

high_b <- list(-0.05,	0.1, 0.1)
med_b <- list(-0.475, 1, 1)
low_b <- list(-1.43, 3, 3)

high_a <- list(-1, 2, 3)
low_a <- list(-0.1, 2, 3)

low_low_g <- list(1.4, 2, 4)
low_high_g <- list(1.2, 2, 4)
high_low_g <- list(3.55, 2, 4)
high_high_g <- list(3.71, 2, 4)

low_p_cens <- 0.25
high_p_cens <- 0.5
