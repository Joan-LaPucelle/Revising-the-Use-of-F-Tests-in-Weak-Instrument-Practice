library(tidyverse)
library(dplyr)
library(mvtnorm)

rm(list=ls())
set.seed(5061)

# Parameters
runs <- 1e6  # number of simulation runs

# Efficient version of table_row_gen
table_row_gen <- function(p, l, rho){
  K <- p
  lambda <- rep(sqrt(l), K)
  constant <- sqrt((1/(rho^2)) - 1)
  
  # Generate all random draws once
  z1 <- mvtnorm::rmvnorm(runs, mean = rep(0, K), sigma = diag(K))
  z2 <- mvtnorm::rmvnorm(runs, mean = rep(0, K), sigma = diag(K))
  nabla <- mvtnorm::rmvnorm(runs, mean = rep(0, K), sigma = diag(K))   # renamed from W
  
  # Precompute denominators: (lambda + z)ᵀ(lambda + z)
  shifted1 <- sweep(z1, 2, lambda, "+")
  denom1   <- rowSums(shifted1^2)
  
  # Relative bias
  num_rel_bias <- rowSums(shifted1 * z1)
  rel_bias_data <- num_rel_bias / denom1
  relative_bias_estimate <- mean(rel_bias_data)
  
  # Second set for absolute error & MSE
  shifted2 <- sweep(z2, 2, lambda, "+")
  denom2   <- rowSums(shifted2^2)
  
  num_bias2 <- rowSums(shifted2 * z2) / denom2
  num_var2  <- rowSums(shifted2 * nabla) / denom2    # now uses nabla
  dataset   <- num_bias2 + constant * num_var2
  
  # Metrics
  relative_absolute_error_estimate <- mean(abs(dataset))
  square_root_of_relative_mse_estimate <- sqrt(mean(dataset^2))
  prob_2SLS_further_than_OLS <- mean(abs(dataset) > 1)
  
  # Output
  vec <- c(rho,
           p,
           l,
           relative_bias_estimate,
           relative_absolute_error_estimate,
           square_root_of_relative_mse_estimate,
           prob_2SLS_further_than_OLS)
  return(vec)
}

# Define parameter vectors
p_vec <- c(4:15, 20, 25, 30, 100, 250)
rho_vec <- rep(0.1, length(p_vec))
l_vec <- c(4.981398, 5.778380, 6.304604, 6.711117, 6.980542, 7.214311,
           7.415126, 7.527178, 7.651759, 7.754555, 7.845581, 7.924440,
           8.185293, 8.338534, 8.448511, 8.834511, 8.933106)

# Run simulations
table <- t(sapply(1:length(p_vec), function(i) {
  table_row_gen(p_vec[i], l_vec[i], rho_vec[i])
}))

# Convert to data frame and rename
table <- as.data.frame(table) %>%
  rename(rho = V1,
         p = V2,
         l_used = V3,
         Relative_Bias_Estimate = V4,
         Relative_Absolute_Error_Estimate = V5,
         Square_Root_of_Relative_Mean_Square_Error_Estimate = V6,
         Prob_2SLS_farther_from_beta0_than_OLS_Estimate = V7)

# Save
saveRDS(table, "Table_for_Paper.rds")
write.csv(table, "Table_for_Paper.csv", row.names = FALSE)