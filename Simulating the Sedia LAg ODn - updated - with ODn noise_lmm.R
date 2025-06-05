simulate_ODn_decay <- function(coef_estimates, coef_se,
                               n_individuals,
                               baseline_mean, baseline_sd,
                               sigma_0, slope_sigma,
                               baseline_noise, fraction,
                               max_follow_up = 10,
                               time_interval = 0.5,
                               noise_model = 1,
                               failure_prob = 0.2,
                               dropout_prob = 0.1) {
  library(truncnorm)
  library(dplyr)
  
  full_time_points <- seq(0, max_follow_up, by = time_interval)
  
  baselines <- truncnorm::rtruncnorm(n_individuals, mean = baseline_mean, sd = baseline_sd, a = 0.03, b = 7.4)
  
  fail_flags <- runif(n_individuals) < failure_prob
  fail_times <- ifelse(fail_flags, runif(n_individuals, min = 1, max = max_follow_up), NA)
  
  dropout_flags <- ifelse(!fail_flags, runif(n_individuals) < dropout_prob, FALSE)
  dropout_times <- ifelse(dropout_flags, runif(n_individuals, min = 1, max = max_follow_up), NA)
  
  # Fixed decay rate per individual
  generate_decay_rate <- function() {
    a <- rnorm(1, mean = coef_estimates[[1]], sd = coef_se[[1]])
    b <- rnorm(1, mean = coef_estimates[[2]], sd = coef_se[[2]])
    decay_rate <- max(a + b * 0, 0.05)  # fixed at baseline
    return(c(decay_rate, a, b))
  }
  
  # Nonlinear decay function
  exp_decay <- function(t, baseline, decay_rate) {
    pmax(exp(intercept + slope * t), 0.001)#pmax(baseline / (1 + decay_rate * t), 0.001)
  }
  
  exp_rebound <- function(start_value, t, rebound_rate) {
    rebound <- start_value * exp(rebound_rate * (t - t[1]))
    pmax(rebound, 0.001)
  }
  
  compute_noise_sd1 <- function(odn) pmax(sigma_0 + slope_sigma * odn, 0.01)
  compute_noise_sd2 <- function(odn) pmax(baseline_noise + fraction * odn, 0.01)
  
  decay_data_list <- vector("list", n_individuals)
  decay_params <- matrix(NA, n_individuals, 2)
  rebound_rate_vector <- rep(NA, n_individuals)
  censor_time_vec <- rep(NA, n_individuals)
  treatment_failure_vector <- rep(FALSE, n_individuals)
  follow_up_years <- rep(NA, n_individuals)
  
  for (i in 1:n_individuals) {
    baseline <- baselines[i]
    has_failed <- fail_flags[i]
    fail_time <- fail_times[i]
    
    if (has_failed) {
      times <- seq(0, fail_time + time_interval, by = time_interval)
      censor_time <- max(times)
    } else if (dropout_flags[i]) {
      times <- seq(0, dropout_times[i], by = time_interval)
      censor_time <- dropout_times[i]
    } else {
      times <- full_time_points
      censor_time <- max_follow_up
    }
    
    follow_up_years[i] <- max(times)
    
    decay_out <- generate_decay_rate()
    decay_rate <- decay_out[1]
    a <- decay_out[2]
    b <- decay_out[3]
    decay_params[i, ] <- c(a, b)
    
    if (has_failed) {
      fail_idx <- which(times >= fail_time)
      treatment_failure_vector[i] <- TRUE
      rebound_rate <- runif(1, min = 0.01, max = 0.05)
      rebound_rate_vector[i] <- rebound_rate
      
      pre_fail_times <- times[1:(fail_idx[1] - 1)]
      expected_odn <- exp_decay(pre_fail_times, baseline, decay_rate)
      
      rebound_t <- times[fail_idx[1]:(fail_idx[1])]
      odn_start <- tail(expected_odn, 1)
      rebound_odn <- exp_rebound(odn_start, rebound_t, rebound_rate)
      
      combined_odn <- c(expected_odn, rebound_odn)
    } else {
      combined_odn <- exp_decay(times, baseline, decay_rate)
    }
    
    noise_sd <- if (noise_model == 1) compute_noise_sd1(combined_odn) else compute_noise_sd2(combined_odn)
    noisy_odn <- pmax(rnorm(length(times), mean = combined_odn, sd = noise_sd), 0.001)
    
    decay_data_list[[i]] <- data.frame(
      individual = i,
      time = times,
      value = noisy_odn
    )
    
    censor_time_vec[i] <- censor_time
  }
  
  decay_data <- do.call(rbind, decay_data_list)
  
  model_parameters <- data.frame(
    id = 1:n_individuals,
    baseline = baselines,
    a = decay_params[, 1],
    b = decay_params[, 2],
    follow_up_years = follow_up_years,
    dropout_time = dropout_times,
    dropped_out = dropout_flags,
    treatment_failure = treatment_failure_vector,
    failure_time = fail_times,
    rebound_rate = rebound_rate_vector,
    censor_time = censor_time_vec
  )
  
  return(list(
    decay_data = decay_data,
    model_parameters = model_parameters
  ))
}

# using actual model formulation
simulate_ODn_decay <- function(coef_estimates, coef_se,
                               n_individuals,
                               baseline_mean, baseline_sd,
                               sigma_0, slope_sigma,
                               baseline_noise, fraction,
                               max_follow_up = 10,
                               time_interval = 0.5,
                               noise_model = 1,
                               failure_prob = 0.2,
                               dropout_prob = 0.1) {
  library(truncnorm)
  library(dplyr)
  
  full_time_points <- seq(0, max_follow_up, by = time_interval)
  
  # Truncated normal baseline ODn values
  baselines <- truncnorm::rtruncnorm(n_individuals, mean = baseline_mean, sd = baseline_sd, a = 0.03, b = 7.4)
  
  # Treatment failure and dropout indicators
  fail_flags <- runif(n_individuals) < failure_prob
  fail_times <- ifelse(fail_flags, runif(n_individuals, min = 1, max = max_follow_up), NA)
  dropout_flags <- ifelse(!fail_flags, runif(n_individuals) < dropout_prob, FALSE)
  dropout_times <- ifelse(dropout_flags, runif(n_individuals, min = 1, max = max_follow_up), NA)
  
  # Decay function (exponential decay)
  exp_decay <- function(t, baseline, decay_rate) {
    pmax(baseline * exp(-decay_rate * t), 0.001)
  }
  
  # Noise models
  compute_noise_sd1 <- function(odn) pmax(sigma_0 + slope_sigma * odn, 0.01)
  compute_noise_sd2 <- function(odn) pmax(baseline_noise + fraction * odn, 0.01)
  
  # Output containers
  decay_data_list <- vector("list", n_individuals)
  decay_params <- matrix(NA, n_individuals, 2)
  censor_time_vec <- rep(NA, n_individuals)
  treatment_failure_vector <- rep(FALSE, n_individuals)
  rebound_rate_vector <- rep(NA, n_individuals)
  follow_up_years <- rep(NA, n_individuals)
  
  for (i in 1:n_individuals) {
    baseline <- baselines[i]
    has_failed <- fail_flags[i]
    fail_time <- fail_times[i]
    
    # Determine time points
    if (has_failed) {
      times <- seq(0, fail_time + time_interval, by = time_interval)
      censor_time <- max(times)
    } else if (dropout_flags[i]) {
      times <- seq(0, dropout_times[i], by = time_interval)
      censor_time <- dropout_times[i]
    } else {
      times <- full_time_points
      censor_time <- max_follow_up
    }
    
    follow_up_years[i] <- max(times)
    
    # Draw log-linear decay parameters and exponentiate for exponential model
    a_log <- rnorm(1, mean = coef_estimates[[1]], sd = coef_se[[1]])
    b_log <- rnorm(1, mean = coef_estimates[[2]], sd = coef_se[[2]])
    decay_rate <- exp(a_log + b_log * 0)  # At baseline
    
    decay_params[i, ] <- c(a_log, b_log)
    
    # Simulate ODn values
    expected_odn <- exp_decay(times, baseline, decay_rate)
    
    # Add noise
    noise_sd <- if (noise_model == 1) compute_noise_sd1(expected_odn) else compute_noise_sd2(expected_odn)
    noisy_odn <- pmax(rnorm(length(times), mean = expected_odn, sd = noise_sd), 0.001)
    
    # Save data
    decay_data_list[[i]] <- data.frame(
      individual = i,
      time = times,
      value = noisy_odn
    )
    
    treatment_failure_vector[i] <- has_failed
    censor_time_vec[i] <- censor_time
  }
  
  decay_data <- do.call(rbind, decay_data_list)
  
  model_parameters <- data.frame(
    id = 1:n_individuals,
    baseline = baselines,
    a = decay_params[, 1],
    b = decay_params[, 2],
    follow_up_years = follow_up_years,
    dropout_time = dropout_times,
    dropped_out = dropout_flags,
    treatment_failure = treatment_failure_vector,
    failure_time = fail_times,
    rebound_rate = rebound_rate_vector,
    censor_time = censor_time_vec
  )
  
  return(list(
    decay_data = decay_data,
    model_parameters = model_parameters
  ))
}


result <- simulate_ODn_decay(
  coef_estimates = coef_estimates,
  coef_se = coef_se,
  n_individuals = 10000,
  baseline_mean = 3.47,
  baseline_sd = 1.55,
  sigma_0 = -0.01469,
  slope_sigma = 0.14513,
  baseline_noise = 0.1,
  fraction = 0.15,
  noise_model = 2  # 1 = heteroskedastic, 2 = realistic noise
)

decay_data <- result$decay_data
model_parameters <- result$model_parameters


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

saveRDS(model_parameters, 'data/Exponential.rds')
jpeg('other_figures/simulated_plot - exponential_new.jpeg', units = "in", width = 9, height = 9, res = 300)
ggplot_plots
dev.off()

# # Load necessary libraries
# library(ggplot2)
# library(dplyr)
# 
# # --- Simulate Data ---
# set.seed(2025)
# 
# sim <- simulate_ODn_decay(
#   coef_estimates = c(0.5, 0.1),   # example values
#   coef_se = c(0.05, 0.02),
#   n_individuals = 100,
#   baseline_mean = 3,
#   baseline_sd = 1,
#   sigma_0 = 0.2,
#   slope_sigma = 0.1,
#   baseline_noise = 0.2,
#   fraction = 0.1,
#   max_follow_up = 10,
#   time_interval = 0.5,
#   noise_model = 1,
#   failure_prob = 0.2,
#   dropout_prob = 0.1
# )
# 
# decay_data <- sim$decay_data
# model_parameters <- sim$model_parameters
# 
# # --- Plot 1: All Individual Trajectories ---
# ggplot(decay_data, aes(x = time, y = value, group = individual)) +
#   geom_line(alpha = 0.3) +
#   labs(
#     title = "ODn Trajectories Over Time",
#     subtitle = "All individuals (n = 100)",
#     x = "Time since ART start (years)",
#     y = "ODn Value"
#   ) +
#   theme_minimal()
# 
# # --- Plot 2: Highlight Failures vs. Non-failures ---
# decay_data <- decay_data %>%
#   left_join(model_parameters %>% select(id, treatment_failure), by = c("individual" = "id"))
# 
# ggplot(decay_data, aes(x = time, y = value, group = individual, color = treatment_failure)) +
#   geom_line(alpha = 0.6) +
#   scale_color_manual(values = c("gray40", "red")) +
#   labs(
#     title = "ODn Trajectories with Treatment Failures Highlighted",
#     x = "Time since ART start (years)",
#     y = "ODn Value",
#     color = "Treatment Failure"
#   ) +
#   theme_minimal()
# 
# # --- Plot 3: Smoothed Trend ---
# ggplot(decay_data, aes(x = time, y = value)) +
#   geom_smooth(method = "loess", span = 0.3, se = FALSE, color = "blue") +
#   labs(
#     title = "Smoothed Average ODn Decay Over Time",
#     x = "Time since ART start (years)",
#     y = "Mean ODn Value"
#   ) +
#   theme_minimal()
# 
# # --- Plot 4: Histogram of Follow-up Duration ---
# ggplot(model_parameters, aes(x = follow_up_years)) +
#   geom_histogram(binwidth = 0.5, fill = "steelblue", color = "white") +
#   labs(
#     title = "Distribution of Follow-up Duration",
#     x = "Follow-up Time (years)",
#     y = "Number of Individuals"
#   ) +
#   theme_minimal()

simulate_ODn_decay_LMM <- function(n_individuals,
                                   baseline_mean, baseline_sd,
                                   slope_mean, slope_sd,
                                   sigma_0 = 0.1, slope_sigma = 0.05,
                                   baseline_noise = 0.05, fraction = 0.1,
                                   max_follow_up = 10,
                                   time_interval = 0.5,
                                   noise_model = 1,
                                   failure_prob = 0.2,
                                   dropout_prob = 0.1) {
  library(truncnorm)
  library(dplyr)
  set.seed(11)
  full_time_points <- seq(0, max_follow_up, by = time_interval)
  
  # Baseline ODn values (truncated normal to avoid negative values)
  baselines <- truncnorm::rtruncnorm(n_individuals, mean = baseline_mean, sd = baseline_sd, a = 0.1, b = 7.4)
  
  # Random individual slopes (may be positive or negative)
  slopes <- rnorm(n_individuals, mean = slope_mean, sd = slope_sd)
  
  # Treatment failure and dropout flags
  fail_flags <- runif(n_individuals) < failure_prob
  fail_times <- ifelse(fail_flags, runif(n_individuals, min = 1, max = max_follow_up), NA)
  
  dropout_flags <- ifelse(!fail_flags, runif(n_individuals) < dropout_prob, FALSE)
  dropout_times <- ifelse(dropout_flags, runif(n_individuals, min = 1, max = max_follow_up), NA)
  
  compute_noise_sd1 <- function(odn) pmax(sigma_0 + slope_sigma * odn, 0.01)
  compute_noise_sd2 <- function(odn) pmax(baseline_noise + fraction * odn, 0.01)
  
  decay_data_list <- vector("list", n_individuals)
  rebound_rate_vector <- rep(NA, n_individuals)
  follow_up_years <- rep(NA, n_individuals)
  censor_time_vec <- rep(NA, n_individuals)
  
  for (i in 1:n_individuals) {
    baseline <- baselines[i]
    slope <- slopes[i]
    
    if (fail_flags[i]) {
      times <- seq(0, fail_times[i] + time_interval, by = time_interval)
      censor_time <- max(times)
    } else if (dropout_flags[i]) {
      times <- seq(0, dropout_times[i], by = time_interval)
      censor_time <- dropout_times[i]
    } else {
      times <- full_time_points
      censor_time <- max_follow_up
    }
    
    follow_up_years[i] <- max(times)
    
    if (fail_flags[i]) {
      fail_idx <- which(times >= fail_times[i])
      rebound_rate <- runif(1, min = 0.01, max = 0.05)
      rebound_rate_vector[i] <- rebound_rate
      
      pre_fail_times <- times[1:(fail_idx[1] - 1)]
      post_fail_times <- times[fail_idx[1]:length(times)]
      
      expected_pre_fail <- baseline + slope * pre_fail_times
      start_value <- tail(expected_pre_fail, 1)
      expected_post_fail <- start_value * exp(rebound_rate * (post_fail_times - post_fail_times[1]))
      
      combined_odn <- c(expected_pre_fail, expected_post_fail)
    } else {
      combined_odn <- baseline + slope * times
    }
    
    combined_odn <- pmax(combined_odn, 0.001)
    
    noise_sd <- if (noise_model == 1) compute_noise_sd1(combined_odn) else compute_noise_sd2(combined_odn)
    noisy_odn <- pmax(rnorm(length(times), mean = combined_odn, sd = noise_sd), 0.001)
    
    decay_data_list[[i]] <- data.frame(
      individual = i,
      time = times,
      value = noisy_odn
    )
    
    censor_time_vec[i] <- censor_time
  }
  
  decay_data <- do.call(rbind, decay_data_list)
  
  model_parameters <- data.frame(
    id = 1:n_individuals,
    baseline = baselines,
    slope = slopes,
    follow_up_years = follow_up_years,
    dropout_time = dropout_times,
    dropped_out = dropout_flags,
    treatment_failure = fail_flags,
    failure_time = fail_times,
    rebound_rate = rebound_rate_vector,
    censor_time = censor_time_vec
  )
  
  return(list(
    decay_data = decay_data,
    model_parameters = model_parameters
  ))
}

sim_results <- simulate_ODn_decay_LMM(
  n_individuals = 1000,
  baseline_mean = 3.0,
  baseline_sd = 0.8,
  slope_mean = -0.083, #-0.25,
  slope_sd = 0.1,
  # sigma_noise = 0.3,        # fixed Gaussian noise SD
  failure_prob = 0.2,
  dropout_prob = 0.1
)

decay_data <-sim_results$decay_data
model_parameters <- sim_results$model_parameters

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

saveRDS(model_parameters, 'data/Exponential_lmm.rds')
ggplot_plots


best_model_choice_LMM <- function(test_data, param_sim_data) {
  time_vec <- test_data[, 2]
  ODn_vec <- test_data[, 3]
  
  dt <- param_sim_data
  baseline <- dt$baseline
  slope <- dt$slope
  
  # browser()
  # Loop over each timepoint
  for (i in 1:nrow(test_data)) {
    t <- time_vec[[i, 1]]
    y <- ODn_vec[[i, 1]]
    pred <- pmax((baseline + slope * t), 0.001)
    dt[[paste0("squared_error", i)]] <- (y - pred)^2
    dt[[paste0("abs_error", i)]] <- abs(y - pred)
  }
  
  # Filter out rows with non-finite values
  dt_clean <- dt[!is.infinite(rowSums(dt %>% dplyr::select(starts_with("squared_error")))), ]
  # browser()
  # Compute RMSE and MAE
  dt_scores <- dt_clean %>%
    rowwise() %>%
    mutate(
      rmse = sqrt(mean(c_across(starts_with("squared_error")))),
      mae = mean(c_across(starts_with("abs_error"))),
      complexity = abs(slope) #+ abs(c)
    ) %>%
    ungroup()
  
  # Step 1: Filter for lowest RMSE
  min_rmse <- min(dt_scores$rmse)
  dt_min_rmse <- dt_scores %>% filter(rmse == min_rmse)
  
  # Step 2: If ties, use MAE
  if (nrow(dt_min_rmse) > 1) {
    min_mae <- min(dt_min_rmse$mae)
    dt_min_rmse <- dt_min_rmse %>% filter(mae == min_mae)
  }
  
  # Step 3: If still ties, use lowest model complexity
  if (nrow(dt_min_rmse) > 1) {
    min_complexity <- min(dt_min_rmse$complexity)
    dt_min_rmse <- dt_min_rmse %>% filter(complexity == min_complexity)
  }
  
  # Step 4: Still ties? Pick one randomly
  if (nrow(dt_min_rmse) > 1) {
    set.seed(11)
    dt_min_rmse <- dt_min_rmse[sample(1:nrow(dt_min_rmse), 1), ]
  }
  
  # Final selection
  return(dt_min_rmse %>% dplyr::select(id, baseline, slope, rmse, mae))
}


results <- list()
ids <- unique(test_data$record_id)

for (i in seq_along(ids)) {
  individual_data <- test_data %>%
    filter(is.na(peak_visit)) %>%
    filter(record_id == ids[i])
  
  best_model_parameters <- best_model_choice_LMM(
    test_data = individual_data,
    param_sim_data = model_parameters#,
    # baselines = model_parameters$baselines
  )
  
  results[[i]] <- best_model_parameters
}

results_df <- do.call(rbind, results)

results <- results_df %>%
  mutate(record_id = ids)

compare_value_with_others_LMM <- function(data_set, t, y_ODn, sigma_ODn, sigma_y_ODn, id) {
  t <- t
  y <- y_ODn
  
  # Compute predicted ODn value
  y_hat <- with(data_set, pmax(baseline + slope * t, 0.001))
  
  # Compute pooled standard deviation
  sigma_pooled <- sqrt(sigma_y_ODn^2 + sigma_ODn^2)
  
  # Compute z-statistic
  z_stat <- (y - y_hat) / sigma_pooled
  
  # One-sided p-value for H1: y > y_hat
  p_value <- 1 - pnorm(z_stat)
  
  # Flag values where y > y_hat is statistically significant at alpha = 0.05
  flagged <- p_value < 0.05
  
  # Return a data frame with full diagnostics
  results <- data.frame(
    id = data_set$id,
    value = y,
    predicted = y_hat,
    z_stat = z_stat,
    p_value = p_value,
    flagged = flagged,
    record_id = id
  )
  
  return(results)
}

test_data_last_visit <- test_data %>%
  filter(!is.na(peak_visit))
results1 <- c()
for (i in 1:length(results$record_id)) {
  results1 <- rbind(results1, compare_value_with_others_LMM(data_set = results %>% filter(record_id == results$record_id[i]),#best_model_parameters, 
                                                        t = test_data_last_visit$time_vec[test_data_last_visit$record_id == results$record_id[i]], #Must be in days
                                                        y_ODn = test_data_last_visit$ODn_vec[test_data_last_visit$record_id == results$record_id[i]], 
                                                        sigma_ODn = sd(test_data$ODn_vec[test_data$record_id == results$record_id[i] & is.na(test_data$peak_visit)]),
                                                        sigma_y_ODn = noise$coefficients[1,1] + noise$coefficients[2,1] * (test_data_last_visit$time_vec[test_data_last_visit$record_id == results$record_id[i]]), #sd((full_dataset %>% filter(Group == 'early suppressions'))$sedia_ODn), #sd(dt$ODn),
                                                        id = results$record_id[i]
  ))
}

dt05 <- as.data.frame(results1) %>%
  left_join(bind_rows(
    dt02 %>%
      mutate(subject_label_blinded = as.character(subject_label_blinded)), dt04
  ) %>%
    mutate(record_id = as.numeric(as.factor(subject_label_blinded))), by = 'record_id') %>%
  dplyr::select(names(as.data.frame(results1)), strata) %>%
  distinct(record_id, .keep_all = T) %>%
  mutate(y_hat_status = as.factor(ifelse(p_value < 0.05, 1, 2)),
         y = as.factor(ifelse(strata == 'early suppressions', 1, 2)) )

dt06 <- dt05 %>%
  dplyr::select(strata, y_hat_status) %>%
  tbl_summary(by = strata
  )
dt06