
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

