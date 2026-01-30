library(tidyverse)
library(dplyr)
library(mvtnorm)

rm(list = ls())

set.seed(2609)

# Inputs ----------
p_values <- c(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3)  # descending order
rho_vec <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
runs <- 1e3

# Function that computes critical value for given (p, l)
func <- function(p, l, z_v, nabla) {
  lambda <- rep(sqrt(l), p)
  lambda_mat <- matrix(lambda, nrow = runs, ncol = p, byrow = TRUE)
  
  ZplusL <- lambda_mat + z_v
  base <- rowSums(ZplusL^2)
  
  prob_est_minus_5_percent <- function(critical_value) {
    mean(base > (critical_value * p)) - 0.05
  }
  
  result <- uniroot(prob_est_minus_5_percent, lower = 0, upper = 100)
  return(result$root)
}

# Compute results_df (first stage, as in your original code)
l_outputs_for_p <- function(p){
  start_time <- Sys.time()
  cat("Starting l_outputs_for_p for p =", p, "\n")
  
  z_v <- mvtnorm::rmvnorm(runs, mean = rep(0,p), sigma = diag(p))
  nabla <- mvtnorm::rmvnorm(runs, mean = rep(0,p), sigma = diag(p))
  
  relative_mean_square_error <- function(l, rho){
    lambda <- rep(sqrt(l), p)
    constant <- (1/(rho^2)) - 1
    lambda_mat <- matrix(lambda, nrow = runs, ncol = p, byrow = TRUE)
    ZplusL <- lambda_mat + z_v
    
    base <- rowSums(ZplusL^2)
    bias_numerator <- rowSums(ZplusL * z_v)
    var_numerator  <- rowSums(ZplusL * nabla)
    
    rel_mean_square_error_data <- (bias_numerator/base)^2 + 
      constant * (var_numerator/base)^2
    
    mean(rel_mean_square_error_data)
  }
  
  objective_function <- function(l, rho) {
    abs(relative_mean_square_error(l, rho) - 1)
  }
  
  l_search <- function(rho_value){
    optim_result <- optim(
      par = 1, fn = objective_function,
      rho = rho_value,
      method = "Brent", lower = 0, upper = 50
    )
    return(as.numeric(optim_result$par))
  }
  
  res <- sapply(rho_vec, l_search)
  end_time <- Sys.time()
  cat("Finished l_outputs_for_p for p =", p, "in", 
      round(difftime(end_time, start_time, units="mins"),2), "minutes\n")
  return(res)
}

# Stage 1: build results_df
cat("Stage 1: Building results_df...\n")
stage1_start <- Sys.time()
results_list <- lapply(p_values, l_outputs_for_p)
results_mat <- do.call(rbind, results_list)
results_df <- data.frame(p = p_values, results_mat)
colnames(results_df) <- c("p", paste0("Rho equals ", rho_vec))
stage1_end <- Sys.time()
cat("Stage 1 finished in", round(difftime(stage1_end, stage1_start, units="mins"),2), "minutes\n")

# Stage 2: transform results_df with func
cat("Stage 2: Transforming results_df...\n")
stage2_start <- Sys.time()
results_df_transformed <- results_df

for (i in seq_len(nrow(results_df_transformed))) {
  p_val <- results_df_transformed$p[i]
  cat("Row", i, "/", nrow(results_df_transformed), " (p =", p_val, ") starting...\n")
  row_start <- Sys.time()
  
  z_v <- mvtnorm::rmvnorm(runs, mean = rep(0,p_val), sigma = diag(p_val))
  nabla <- mvtnorm::rmvnorm(runs, mean = rep(0,p_val), sigma = diag(p_val))
  
  for (j in 2:ncol(results_df_transformed)) {
    l_val <- results_df_transformed[i, j]
    results_df_transformed[i, j] <- func(p_val, l_val, z_v, nabla)
    cat("  Finished col", j-1, "/", ncol(results_df_transformed)-1, "for p =", p_val, "\n")
  }
  
  row_end <- Sys.time()
  cat("Row", i, "completed in", 
      round(difftime(row_end, row_start, units="mins"),2), "minutes\n")
}

stage2_end <- Sys.time()
cat("Stage 2 finished in", round(difftime(stage2_end, stage2_start, units="mins"),2), "minutes\n")

print(results_df_transformed)

# Invert the row order so that p=3 comes first and p=15 comes last
results_df_transformed_inverted <- results_df_transformed %>%
  arrange(p)

# Save outputs
saveRDS(results_df_transformed_inverted, "FinalTableCriticalValues.rds")
write.csv(results_df_transformed_inverted, "FinalTableCriticalValues.csv", row.names = FALSE)
