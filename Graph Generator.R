library(tidyverse)
library(dplyr)
library(mvtnorm)
library(ggsci)
library(reshape2)

rm(list = ls())
set.seed(2609)

# Parameters
runs <- 1000000
p <- 5

# Generate Monte Carlo Samples once
z_v <- mvtnorm::rmvnorm(runs, mean = rep(0, p), sigma = diag(p))
w   <- mvtnorm::rmvnorm(runs, mean = rep(0, p), sigma = diag(p))

# Vectorized raw_data_gen
raw_data_gen <- function(l, rho){
  lambda   <- rep(sqrt(l), p)
  constant <- (1 / (rho^2)) - 1
  
  shifted <- sweep(z_v, 2, lambda, "+")   # z + lambda
  
  base           <- rowSums(shifted^2)
  bias_numerator <- rowSums(shifted * z_v)
  var_numerator  <- rowSums(shifted * w)
  
  # Probability
  prob_est <- as.numeric(base > (10 * p))
  prob_F_greater_than_ten_est <- mean(prob_est)
  prob_F_greater_than_ten_est_std_error <- sd(prob_est)
  
  # Relative bias
  rel_bias_data <- bias_numerator / base
  rel_bias_est <- mean(rel_bias_data)
  rel_bias_est_std_error <- sd(rel_bias_data)
  
  # Relative absolute error
  rel_abs_error_data <- abs((bias_numerator + sqrt(constant) * var_numerator) / base)
  rel_abs_error_est <- mean(rel_abs_error_data)
  rel_abs_error_est_std_error <- sd(rel_abs_error_data)
  
  # Relative mean square error
  rel_mean_square_error_data <- (bias_numerator/base)^2 + constant * ((var_numerator/base)^2)
  rel_mean_square_error_est <- mean(rel_mean_square_error_data)
  rel_mean_square_error_est_std_error <- sd(rel_mean_square_error_data)
  sqrt_rel_mean_square_error_est <- sqrt(rel_mean_square_error_est)
  
  vec <- c(l,
           prob_F_greater_than_ten_est,
           prob_F_greater_than_ten_est_std_error,
           rel_bias_est,
           rel_bias_est_std_error,
           rel_abs_error_est, 
           rel_abs_error_est_std_error,
           sqrt_rel_mean_square_error_est,
           rel_mean_square_error_est_std_error, 
           rho)
  
  return(as.matrix(vec))
}

# Input vectors
rho_vec <- rep(0.1, 126)
rho_vec[22:42]  <- 0.3
rho_vec[43:63]  <- 0.5
rho_vec[64:84]  <- 0.7
rho_vec[85:105] <- 0.9
rho_vec[106:126]<- 1.0

l_options <- seq(0, 10, by = 0.5)
l_vec <- rep(l_options, times = 6)
length <- length(l_vec)

# Run simulations
table <- sapply(1:length, function(s) {
  raw_data_gen(l_vec[s], rho_vec[s])
})

table <- t(table) |> as.data.frame()

# Rename columns
table <- table %>% rename(
  l = V1,
  Prob_F_greater_than_10_estimate = V2,
  Prob_F_greater_than_10_Estimate_std_error = V3,
  Relative_Bias_Estimate = V4, 
  Relative_Bias_Estimate_std_error = V5,
  Relative_Absolute_Error_Estimate = V6, 
  Relative_Absolute_Error_Estimate_std_error = V7,
  Square_Root_of_Relative_Mean_Square_Error_Estimate = V8, 
  Relative_Mean_Square_Error_Estimate_std_error = V9, 
  Rho = V10
)

saveRDS(table, "Table for Graphs.rds")
write.csv(table, "Table for Graphs.csv")

#### Plot Graphs ---------

dataset <- table 
dataset <- dataset %>% select(
  -Relative_Mean_Square_Error_Estimate_std_error, 
  -Relative_Absolute_Error_Estimate_std_error, 
  -Relative_Bias_Estimate_std_error, 
  -Prob_F_greater_than_10_Estimate_std_error
)

dataset <- melt(dataset, na.rm = FALSE, value.name = "Metric_Value",
                id = c("l", "Rho", "Prob_F_greater_than_10_estimate")) %>% 
  mutate(Metric = case_when(
    variable == "Relative_Bias_Estimate" ~ "Relative Bias", 
    variable == "Relative_Absolute_Error_Estimate" ~ "Relative Absolute Error", 
    variable == "Square_Root_of_Relative_Mean_Square_Error_Estimate" ~ "Square Root of RMSE"
  )) %>%
  select(-variable)

### Plot 1: Metrics vs Expected F Value (for fixed rho = 0.1)
data_plot <- dataset %>% filter(Rho == 0.1)

# Add Expected F Value = 1 + l
data_plot <- data_plot %>%
  mutate(Expected_F_Value = 1 + l) %>%
  filter(!(Prob_F_greater_than_10_estimate <= 0.000002 &
             Metric %in% c("Relative Absolute Error", "Square Root of RMSE")))

# Find approximate Expected F Values for Prob(F > 10) ≈ 0.05, 0.25, 0.5
ref_probs <- c(0.05, 0.25, 0.5)
ref_points <- sapply(ref_probs, function(prob) {
  subset_data <- data_plot %>%
    filter(abs(Prob_F_greater_than_10_estimate - prob) == 
             min(abs(Prob_F_greater_than_10_estimate - prob)))
  unique(subset_data$Expected_F_Value)
})

# Plot with vertical lines at reference Expected F Values
p1 <- ggplot(data = data_plot, aes(x = Expected_F_Value,
                                   y = Metric_Value, 
                                   shape = Metric, 
                                   colour = Metric)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = ref_points, linetype = "dotted", linewidth = 0.8, color = "black") +
  annotate("text", x = ref_points, 
           y = max(data_plot$Metric_Value) * 0.8,  # lowered text position
           label = paste0("P(F>10) = ", ref_probs),
           angle = 90, vjust = -0.5, hjust = 1.1, size = 3.5, color = "black") +
  theme_bw(base_size = 14) +
  scale_color_lancet() +
  labs(
    x = "Expected Value of F Statistic",
    y = "Metric Value", 
    title = "Plot Metrics against Expected F Value (ρ = 0.1)"
  ) +
  theme(
    legend.position = c(0.83, 0.9),  # moved further upward & right
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    legend.title = element_blank()
  )

ggsave("metrics_vs_expectedF.png", plot = p1, width = 8, height = 6, dpi = 300)

### Plot 2: Metrics vs Rho (fixed l = 6)
data_plot <- dataset %>% filter(l == 6)

p2 <- ggplot(data = data_plot, aes(x = Rho,
                                   y = Metric_Value, 
                                   shape = Metric, 
                                   colour = Metric)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 14) +
  scale_color_lancet() +
  labs(
    x = "Rho",
    y = "Metric Value", 
    title = "Plot Metrics against Rho"
  ) +
  theme(
    legend.position = c(0.75, 0.75),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    legend.title = element_blank()
  )

ggsave("metrics_vs_rho.png", plot = p2, width = 8, height = 6, dpi = 300)
