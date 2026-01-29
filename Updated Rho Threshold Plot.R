# Efficient Rho Threshold Simulation and Plotting
# X-axis = Expectation of F Statistic (E[F] = 1 + l)
# Dotted vertical lines at min E[F] such that P(F>10) >= {0.05, 0.25, 0.5}
# (minimized over all p shown in each plot)

library(tidyverse)
library(dplyr)
library(ggsci)
library(reshape2)
library(mvtnorm)

rm(list = ls())
set.seed(2609)

# ----------------- INPUTS -----------------
runs <- 1000000   # number of simulation runs
p_vec <- c(5, 6, 7, 8, 9, 10, 15, 20)
l_vec <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5,
           6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12, 13, 14)

# ----------------- SIMULATION FUNCTION -----------------
simulate_for_p <- function(p) {
  message("Running simulations for p = ", p)
  
  # Generate Monte Carlo Samples
  z_v <- rmvnorm(runs, mean = rep(0, p), sigma = diag(p))
  w   <- rmvnorm(runs, mean = rep(0, p), sigma = diag(p))
  
  raw_data_gen <- function(l) {
    lambda <- rep(sqrt(l), p)
    
    shifted <- sweep(z_v, 2, lambda, "+")     # z + lambda
    base <- rowSums(shifted^2)
    
    bias_numerator <- rowSums(shifted * z_v)
    var_numerator  <- rowSums(shifted * w)
    
    # P(F > 10) where F = base/p, so F>10 <=> base > 10p
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

# ----------------- RUN ALL p VALUES -----------------
results <- lapply(p_vec, simulate_for_p)
names(results) <- paste0("p", p_vec)

# ----------------- SAVE RESULTS -----------------
for (p in p_vec) {
  tbl <- results[[paste0("p", p)]]
  saveRDS(tbl, paste0("Table_Rho_Threshold_p", p, ".rds"))
  write.csv(tbl, paste0("Table_Rho_Threshold_p", p, ".csv"), row.names = FALSE)
}

# ----------------- PLOTTING HELPERS -----------------
# Compute reference x-values: for each probability bound q, take the smallest Expected_F
# among all rows in the plotting dataset such that Prob(F>10) >= q.
get_ref_points_min_expectedF <- function(df, prob_bounds = c(0.05, 0.25, 0.5)) {
  # df must contain: Expected_F, Prob_F_greater_than_10_estimate
  sapply(prob_bounds, function(q) {
    eligible <- df %>% filter(Prob_F_greater_than_10_estimate >= q)
    if (nrow(eligible) == 0) return(NA_real_)
    min(eligible$Expected_F, na.rm = TRUE)
  })
}

# ----------------- PLOTTING: LARGE p -----------------
dataset_large <- bind_rows(results$p10, results$p15, results$p20) %>%
  # keep l to compute Expected_F
  select(-Second_Term_Std_Error, -First_Term_Std_Error,
         -Prob_F_greater_than_10_Estimate_std_error) %>%
  filter(Rho_Threshold <= 0.99, Prob_F_greater_than_10_estimate > 8e-7) %>%
  mutate(
    Expected_F = 1 + l,
    p = factor(p, levels = as.character(c(10, 15, 20)))
  )

ref_probs_large  <- c(0.05, 0.25, 0.5)
ref_points_large <- get_ref_points_min_expectedF(dataset_large, ref_probs_large)

p1 <- ggplot(data = dataset_large,
             aes(x = Expected_F, y = Rho_Threshold, shape = p, colour = p)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = ref_points_large,
             linetype = "dotted", linewidth = 0.8, color = "black") +
  annotate("text",
           x = ref_points_large,
           y = max(dataset_large$Rho_Threshold, na.rm = TRUE) * 0.85,
           label = paste0("P(F>10) \u2265 ", ref_probs_large),
           angle = 90, vjust = -0.5, hjust = 1.1, size = 3.5, color = "black") +
  theme_bw(base_size = 14) +
  scale_color_lancet() +
  labs(
    x = "Expectation of F Statistic",
    y = "Rho Threshold",
    title = "Rho Threshold for when 2SLS is better than OLS",
    shape = "Number of Instruments",
    colour = "Number of Instruments"
  ) +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

ggsave("Rho_Threshold_Large_p.png", plot = p1, width = 8, height = 6, dpi = 300)

# ----------------- PLOTTING: SMALL p -----------------
dataset_small <- bind_rows(results$p5, results$p6, results$p7,
                           results$p8, results$p9, results$p10) %>%
  # keep l to compute Expected_F
  select(-Second_Term_Std_Error, -First_Term_Std_Error,
         -Prob_F_greater_than_10_Estimate_std_error) %>%
  filter(Rho_Threshold <= 0.99, Prob_F_greater_than_10_estimate > 8e-7) %>%
  mutate(
    Expected_F = 1 + l,
    p = factor(p, levels = as.character(5:10))
  )

ref_probs_small  <- c(0.05, 0.25, 0.5)
ref_points_small <- get_ref_points_min_expectedF(dataset_small, ref_probs_small)

p2 <- ggplot(data = dataset_small,
             aes(x = Expected_F, y = Rho_Threshold, shape = p, colour = p)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = ref_points_small,
             linetype = "dotted", linewidth = 0.8, color = "black") +
  annotate("text",
           x = ref_points_small,
           y = max(dataset_small$Rho_Threshold, na.rm = TRUE) * 0.85,
           label = paste0("P(F>10) \u2265 ", ref_probs_small),
           angle = 90, vjust = -0.5, hjust = 1.1, size = 3.5, color = "black") +
  theme_bw(base_size = 14) +
  scale_color_lancet() +
  labs(
    x = "Expectation of F Statistic",
    y = "Rho Threshold",
    title = "Rho Threshold for when 2SLS is better than OLS",
    shape = "Number of Instruments",
    colour = "Number of Instruments"
  ) +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

ggsave("Rho_Threshold_Small_p.png", plot = p2, width = 8, height = 6, dpi = 300)
