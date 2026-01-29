# Efficient Rho Threshold Simulation and Plotting
# Legend placed inside graph

library(tidyverse)
library(dplyr)
library(ggsci)
library(reshape2)
library(mvtnorm)

rm(list = ls())
set.seed(2609)

# Inputs ----------
runs <- 1000000   # number of simulation runs
p_vec <- c(5, 6, 7, 8, 9, 10, 15, 20)
l_vec <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5,
           6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12, 13, 14)

# Function to compute results for one p
simulate_for_p <- function(p) {
  message("Running simulations for p = ", p)
  
  # Generate Monte Carlo Samples
  z_v <- rmvnorm(runs, mean = rep(0, p), sigma = diag(p))
  w   <- rmvnorm(runs, mean = rep(0, p), sigma = diag(p))
  
  raw_data_gen <- function(l) {
    lambda <- rep(sqrt(l), p)
    
    shifted <- sweep(z_v, 2, lambda, "+")
    base <- rowSums(shifted^2)
    
    bias_numerator <- rowSums(shifted * z_v)
    var_numerator  <- rowSums(shifted * w)
    
    prob_est <- as.numeric(base > (10 * p))
    prob_F_greater_than_ten_est <- mean(prob_est)
    prob_F_greater_than_ten_est_std_error <- sd(prob_est)
    
    first_term_data  <- (bias_numerator / base)^2
    second_term_data <- (var_numerator / base)^2
    
    first_term_est  <- mean(first_term_data)
    second_term_est <- mean(second_term_data)
    
    first_term_est_std_error  <- sd(first_term_data)
    second_term_est_std_error <- sd(second_term_data)
    
    rho_threshold <- sqrt(second_term_est / (1 - first_term_est + second_term_est))
    
    return(c(l,
             prob_F_greater_than_ten_est,
             prob_F_greater_than_ten_est_std_error,
             rho_threshold,
             first_term_est_std_error,
             second_term_est_std_error))
  }
  
  table <- t(sapply(l_vec, raw_data_gen)) %>% as.data.frame()
  
  table <- table %>% rename(
    l = V1,
    Prob_F_greater_than_10_estimate = V2,
    Prob_F_greater_than_10_Estimate_std_error = V3,
    Rho_Threshold = V4,
    First_Term_Std_Error = V5,
    Second_Term_Std_Error = V6
  )
  
  table$p <- as.character(p)
  return(table)
}

# Run for all p values
results <- lapply(p_vec, simulate_for_p)
names(results) <- paste0("p", p_vec)

# Save all results
for (p in p_vec) {
  tbl <- results[[paste0("p", p)]]
  saveRDS(tbl, paste0("Table_Rho_Threshold_p", p, ".rds"))
  write.csv(tbl, paste0("Table_Rho_Threshold_p", p, ".csv"), row.names = FALSE)
}

# ----------------- PLOTTING -----------------

# Large p datasets
dataset_large <- bind_rows(results$p10, results$p15, results$p20) %>%
  select(-Second_Term_Std_Error, -First_Term_Std_Error,
         -Prob_F_greater_than_10_Estimate_std_error, -l) %>%
  filter(Rho_Threshold <= 0.99, Prob_F_greater_than_10_estimate > 8e-7)

dataset_large$p <- factor(dataset_large$p, levels = as.character(c(10, 15, 20)))

p1 <- ggplot(data = dataset_large, aes(x = Prob_F_greater_than_10_estimate, 
                                       y = Rho_Threshold,
                                       shape = p, colour = p)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 14) +
  scale_color_lancet() +
  labs(
    x = "Probability F > 10",
    y = "Rho Threshold",
    title = "Rho Threshold for when 2SLS is better than OLS",
    shape = "Number of Instruments",
    colour = "Number of Instruments"
  ) +
  theme(
    legend.position = c(0.8, 0.8), # Inside top-right
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

ggsave("Rho_Threshold_Large_p.png", plot = p1, width = 8, height = 6, dpi = 300)

# Small p datasets
dataset_small <- bind_rows(results$p5, results$p6, results$p7,
                           results$p8, results$p9, results$p10) %>%
  select(-Second_Term_Std_Error, -First_Term_Std_Error,
         -Prob_F_greater_than_10_Estimate_std_error, -l) %>%
  filter(Rho_Threshold <= 0.99, Prob_F_greater_than_10_estimate > 8e-7)

dataset_small$p <- factor(dataset_small$p, levels = as.character(5:10))

p2 <- ggplot(data = dataset_small, aes(x = Prob_F_greater_than_10_estimate, 
                                       y = Rho_Threshold,
                                       shape = p, colour = p)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 14) +
  scale_color_lancet() +
  labs(
    x = "Probability F > 10",
    y = "Rho Threshold",
    title = "Rho Threshold for when 2SLS is better than OLS",
    shape = "Number of Instruments",
    colour = "Number of Instruments"
  ) +
  theme(
    legend.position = c(0.8, 0.8), # Inside top-right
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

ggsave("Rho_Threshold_Small_p.png", plot = p2, width = 8, height = 6, dpi = 300)