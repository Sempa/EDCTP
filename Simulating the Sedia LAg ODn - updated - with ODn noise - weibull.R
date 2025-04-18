# Load required package
library(truncnorm)
library(dplyr)

# Set seed for reproducibility
set.seed(123)

# Number of individuals
n_individuals <- 10000

# Time points (6-month intervals over 10 years)
time_points <- seq(0, 10, by = 0.5)

# Define the distributions for baseline
baseline_mean <- 3.47
baseline_sd <- 1.55
baselines <- truncnorm::rtruncnorm(n_individuals, mean = baseline_mean, sd = baseline_sd, a = 0.001, b = 5)

# Function to generate Weibull parameters from a two-degree polynomial
generate_weibull_params <- function(t) {
  # Coefficients from polynomial model (example using model summary)
  a <- rnorm(1, mean = summary(poly_model)$tTable[,1][[3]], 
             sd = summary(poly_model)$tTable[,2][[3]] * sqrt(length(unique(sedia_eddi_data$subject_label_blinded))))
  b <- rnorm(1, mean = summary(poly_model)$tTable[,1][[2]], 
             sd = summary(poly_model)$tTable[,2][[2]] * sqrt(length(unique(sedia_eddi_data$subject_label_blinded))))
  c <- rnorm(1, mean = summary(poly_model)$tTable[,1][[1]], 
             sd = summary(poly_model)$tTable[,2][[1]] * sqrt(length(unique(sedia_eddi_data$subject_label_blinded))))
  
  lambda <- pmax(a * t^2 + b * t + c, 0.01)  # Time-varying scale
  k <- runif(1, min = 0.5, max = 1)          # Individual-level shape parameter
  
  return(c(lambda, rep(k, length(t)), a, b, c))  # Return full vector
}

# Generate Weibull parameters for each individual
weibull_params <- sapply(1:n_individuals, function(i) {
  generate_weibull_params(time_points)
})

# Weibull decay function
weibull_decay <- function(t, baseline, lambda, k) {
  pmax(baseline * exp(-(lambda * t)^k), 0.001)
}

# Define measurement noise model
sigma_0 <- -0.01469  # Intercept
slope_sigma <- 0.14513  # Slope

# Function to compute heteroskedastic noise SD based on ODn
compute_noise_sd <- function(odn) {
  pmax(sigma_0 + slope_sigma * odn, 0.01)
}

# Generate decay data with Weibull decay and noise
decay_data <- data.frame(
  time = rep(time_points, n_individuals),
  individual = rep(1:n_individuals, each = length(time_points)),
  value = unlist(lapply(1:n_individuals, function(i) {
    lambda_vec <- weibull_params[1:length(time_points), i]
    k_vec <- weibull_params[(length(time_points) + 1):(2 * length(time_points)), i][1]  # constant k per person
    expected_odn <- weibull_decay(time_points, baselines[i], lambda_vec, k_vec)
    noise_sd <- sapply(expected_odn, compute_noise_sd)
    noisy_odn <- pmax(rnorm(length(expected_odn), mean = expected_odn, sd = noise_sd), 0.001)
    return(noisy_odn)
  }))
)

# Extract polynomial coefficients and shape for each individual
model_parameters <- bind_cols(
  id = 1:n_individuals,
  baselines = baselines,
  a = weibull_params[2 * length(time_points) + 1, ],
  b = weibull_params[2 * length(time_points) + 2, ],
  c = weibull_params[2 * length(time_points) + 3, ],
  shape = weibull_params[length(time_points) + 1, ]  # shape k is the same across time
)
saveRDS(model_parameters, 'data/Weibull.rds')

x <- full_dataset %>%
  group_by(subject_label_blinded) %>%
  arrange(subject_label_blinded, years_since_tx_start) %>%
  mutate(visits = 1:length(subject_label_blinded),
         max_visits = max(visits)) %>%
  filter(Group == 'early suppressions') %>%
  distinct(subject_label_blinded, .keep_all = T)
x2 = cephia_samples %>%
  distinct(subject_label_blinded, .keep_all = T)

dataset_ggplot <- decay_data %>%
  filter(individual %in% sample(length(unique(decay_data$individual)), 100, replace = F)) 
set.seed(11)
ggplot_plots <- ggpubr::ggarrange( 
  ggplot(full_dataset %>%
           filter(Group == 'early suppressions'), 
         aes(x = years_since_tx_start, y = sedia_ODn, group = subject_label_blinded, color = subject_label_blinded)) +
    geom_line(size = 1.5) +
    # geom_line(alpha = 0.5) +
    labs(#title = "Exponential Decay Curves for Individuals",
      x = "Time since ART start (years)",
      y = "ODn Value") +
    # theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      legend.position = "none"
    ),
  
  
  ggplot(decay_data %>%
           mutate(flag = as.logical(ifelse(individual %in% sample(n_individuals, 50, replace = F), 1, 0))), 
         aes(x = time, y = value, group = individual, color = flag)) +
    geom_line(alpha = 0.5, linewidth = 1.5) +
    gghighlight(flag, use_direct_label = FALSE, unhighlighted_colour = "grey70") +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "black")) +
    labs(#title = "Exponential Decay Curves for Individuals",
      x = "Time since ART start (years)",
      y = "Estimated ODn Value")  +
    scale_y_continuous(limits = c(0, max(dataset_ggplot$value))) +
    scale_x_continuous(limits = c(0, max(dataset_ggplot$time))) +
    # theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      legend.position = "none"
    ),
  labels = c("A", "B"),
  ncol = 1, nrow = 2
)

jpeg('other_figures/simulated_plot - weibull.jpeg', units = "in", width = 9, height = 9, res = 300)
ggplot_plots
dev.off()

