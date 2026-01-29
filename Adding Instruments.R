# Optimized version of your simulation code
library(tidyverse)
library(dplyr)
library(mvtnorm)

rm(list=ls())

#Inputs ----------
runs <- 1000000 # number of simulation runs

raw_data_gen <- function(p, l, seed_number){
  set.seed(seed_number)
  
  # Generate Z and W in one shot
  z_v <- mvtnorm::rmvnorm(runs, mean = rep(0, p), sigma = diag(p))
  w   <- mvtnorm::rmvnorm(runs, mean = rep(0, p), sigma = diag(p))
  
  lambda <- rep(sqrt(l), p)
  
  # Precompute shifted z = lambda + z
  z_shift <- sweep(z_v, 2, lambda, "+")
  
  # Base term = (lambda + z)'(lambda + z)
  base <- rowSums(z_shift^2)
  
  # Bias numerator = (lambda + z)'z
  bias_numerator <- rowSums(z_shift * z_v)
  
  # Variance numerator = (lambda + z)'w
  var_numerator <- rowSums(z_shift * w)
  
  # Prob(F > 10p)
  prob_est <- as.numeric(base > (10 * p))
  prob_F_greater_than_ten_est <- mean(prob_est)
  prob_F_greater_than_ten_est_std_error <- sqrt(var(prob_est))
  
  # Relative bias
  rel_bias_data <- bias_numerator / base
  rel_bias_est <- mean(rel_bias_data)
  rel_bias_est_std_error <- sqrt(var(rel_bias_data))
  
  # First RMSE term
  first_term_data <- (bias_numerator / base)^2
  first_term_est <- mean(first_term_data)
  first_term_est_std_error <- sqrt(var(first_term_data))
  
  # Second RMSE term
  second_term_data <- (var_numerator / base)^2
  second_term_est <- mean(second_term_data)
  second_term_est_std_error <- sqrt(var(second_term_data))
  
  # Return vector
  vec <- c(p, 
           l,
           prob_F_greater_than_ten_est,
           prob_F_greater_than_ten_est_std_error,
           rel_bias_est,
           rel_bias_est_std_error,
           first_term_est,
           first_term_est_std_error,
           second_term_est,
           second_term_est_std_error)
  
  return(as.matrix(vec))
}

##ENTER VALUE INPUTS HERE
p_vec    <- c(5, 10, 15, 20, 5, 6, 7, 8, 9, 10)
l_vec    <- c(7, 7, 7, 7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7)
seed_vec <- c(347, 4149, 2428, 4673, 1356, 4497, 235, 1083, 4126, 309)

n_cases <- length(l_vec)

# Run table generation
table <- sapply(1:n_cases, function(s) {
  raw_data_gen(p_vec[s], l_vec[s], seed_vec[s])
})

table <- t(table) %>% 
  as.data.frame() %>% 
  rename(p = V1, 
         l = V2,
         Prob_F_greater_than_10_estimate = V3,
         Prob_F_greater_than_10_Estimate_std_error = V4,
         Relative_Bias_Estimate = V5,
         Relative_Bias_Estimate_Std_Error = V6,
         First_Term_Estimate = V7,
         First_Term_Estimate_Std_Error = V8,
         Second_Term_Estimate = V9,
         Second_Term_Estimate_Std_Error = V10)

# Save outputs
saveRDS(table, "Preliminary Table Adding Instruments Experiment.rds")
write.csv(table, "Preliminary Table Adding Instruments Experiment.csv", row.names = FALSE)
