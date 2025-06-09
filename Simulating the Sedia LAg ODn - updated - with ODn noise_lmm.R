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
    pmax(baseline / (1 + decay_rate * t), 0.001)
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


result <- simulate_ODn_decay1(
  coef_estimates = coef_estimates,
  coef_se = coef_se,
  n_individuals = 10000,
  baseline_mean = 3.47,
  baseline_sd = 1.55,
  sigma_0 = -0.01469,
  slope_sigma = 0.14513,
  baseline_noise = 0.1,
  fraction = 0.15,
  noise_model = 1  # 1 = heteroskedastic, 2 = realistic noise
)

decay_data <- result$decay_data
model_parameters <- result$model_parameters

# using actual model formulation
simulate_ODn_decay1 <- function(coef_estimates, coef_se,
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
  
  # Simulate log-linear coefficients
  generate_log_decay_coefs <- function() {
    a <- rnorm(1, mean = coef_estimates[[1]], sd = coef_se[[1]])
    b <- rnorm(1, mean = coef_estimates[[2]], sd = coef_se[[2]])
    return(c(a, b))  # log(ODn) = a + b*t
  }
  
  # Simulate ODn from log-linear model
  log_linear_decay <- function(t, a, b) {
    log_odn <- a + b * t
    pmax(exp(log_odn), 0.001)
  }
  
  # Rebound function in ODn space
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
    
    coefs <- generate_log_decay_coefs()
    a <- coefs[1]
    b <- coefs[2]
    decay_params[i, ] <- coefs
    
    if (has_failed) {
      fail_idx <- which(times >= fail_time)
      treatment_failure_vector[i] <- TRUE
      rebound_rate <- runif(1, min = 0.01, max = 0.05)
      rebound_rate_vector[i] <- rebound_rate
      
      pre_fail_times <- times[1:(fail_idx[1] - 1)]
      expected_odn <- log_linear_decay(pre_fail_times, a, b)
      
      rebound_t <- times[fail_idx[1]:(fail_idx[1])]
      odn_start <- tail(expected_odn, 1)
      rebound_odn <- exp_rebound(odn_start, rebound_t, rebound_rate)
      
      combined_odn <- c(expected_odn, rebound_odn)
    } else {
      combined_odn <- log_linear_decay(times, a, b)
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

result <- simulate_ODn_decay1(
  coef_estimates = coef_estimates,
  coef_se = coef_se,
  n_individuals = 10000,
  baseline_mean = 3.47,
  baseline_sd = 1.55,
  sigma_0 = -0.01469,
  slope_sigma = 0.14513,
  baseline_noise = 0.1,
  fraction = 0.15,
  noise_model = 1  # 1 = heteroskedastic, 2 = realistic noise
)

decay_data <- result$decay_data
model_parameters <- result$model_parameters

# Adding a sigmoid or logistic rebound
#' Simulate ODn Decay Trajectories with Treatment Failure and Antibody Rebound
#'
#' Simulates individual-level longitudinal ODn decay data using an exponential decay model
#' with optional noise, treatment failure, dropout, and post-failure antibody rebound.
#'
#' @param coef_estimates A numeric vector of length 2 with the mean coefficients (`a`, `b`)
#'   for the log-linear decay model (`log(ODn) = a + b * time`).
#' @param coef_se A numeric vector of length 2 with the standard errors for `a` and `b`.
#' @param n_individuals Integer. Number of individuals to simulate.
#' @param baseline_mean Mean of the baseline ODn distribution (used if truncated normal sampling is enabled).
#' @param baseline_sd Standard deviation of the baseline ODn distribution.
#' @param sigma_0 Intercept for heteroskedastic noise model 1.
#' @param slope_sigma Slope for heteroskedastic noise model 1.
#' @param baseline_noise Intercept for noise model 2.
#' @param fraction Fraction of the current ODn value used to determine standard deviation in noise model 2.
#' @param max_follow_up Numeric. Maximum follow-up duration in years (default is 10).
#' @param time_interval Numeric. Time interval between visits (in years; default is 0.5).
#' @param noise_model Integer. Noise model to use: `1` = sigma depends on `ODn`, `2` = alternative model.
#' @param failure_prob Probability that an individual will experience treatment failure.
#' @param dropout_prob Probability that a non-failing individual will be censored due to dropout.
#' @param rebound_model Character. Either `"sigmoid"` or `"logistic"` to control the rebound trajectory transition.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{decay_data}{A data frame with columns `individual`, `time`, and `value` (simulated ODn values).}
#'   \item{model_parameters}{A data frame with simulation parameters for each individual, including decay coefficients,
#'   baseline ODn, failure/dropout status, censoring time, rebound rate, etc.}
#' }
#'
#' @details
#' Each individual's ODn trajectory follows a log-linear decay. Individuals may experience:
#' \itemize{
#'   \item Treatment failure at a random time (with probability `failure_prob`), after which their ODn rebounds.
#'   \item Dropout (with probability `dropout_prob`) if they do not fail.
#'   \item Rebound is modeled as an exponential increase, blended into the decay using a smooth sigmoid or logistic transition.
#'   \item Gaussian noise is added, with a choice between two heteroskedastic noise models.
#' }
#'
#' The baseline ODn is sampled from a uniform distribution between 0.03 and 7.4. The decay model parameters (`a`, `b`) are sampled per individual.
#'
#' @importFrom truncnorm rtruncnorm
#' @importFrom dplyr %>%
#'
#' @examples
#' sim <- simulate_ODn_decay2(
#'   coef_estimates = c(2, -0.2),
#'   coef_se = c(0.1, 0.05),
#'   n_individuals = 100,
#'   baseline_mean = 2,
#'   baseline_sd = 0.5,
#'   sigma_0 = 0.05,
#'   slope_sigma = 0.1,
#'   baseline_noise = 0.05,
#'   fraction = 0.2,
#'   max_follow_up = 5,
#'   time_interval = 0.5,
#'   noise_model = 1,
#'   failure_prob = 0.3,
#'   dropout_prob = 0.2,
#'   rebound_model = "sigmoid"
#' )
#'
#' head(sim$decay_data)
#' head(sim$model_parameters)
#'
#' @export
simulate_ODn_decay2 <- function(coef_estimates, coef_se,
                                n_individuals,
                                baseline_mean, baseline_sd,
                                sigma_0, slope_sigma,
                                baseline_noise, fraction,
                                max_follow_up = 10,
                                time_interval = 0.5,
                                noise_model = 1,
                                failure_prob = 0.2,
                                dropout_prob = 0.1,
                                rebound_model = c("sigmoid", "logistic")) {
  library(truncnorm)
  library(dplyr)
  
  rebound_model <- match.arg(rebound_model)
  full_time_points <- seq(0, max_follow_up, by = time_interval)
  
  # baselines <- truncnorm::rtruncnorm(n_individuals, mean = baseline_mean, sd = baseline_sd, a = 0.03, b = 7.4)
  baselines <- runif(n_individuals, min = 0.03, max = 7.4)
  fail_flags <- runif(n_individuals) < failure_prob
  fail_times <- ifelse(fail_flags, runif(n_individuals, min = 1, max = max_follow_up - 3 * time_interval), NA)
  dropout_flags <- ifelse(!fail_flags, runif(n_individuals) < dropout_prob, FALSE)
  dropout_times <- ifelse(dropout_flags, runif(n_individuals, min = 1, max = max_follow_up), NA)
  
  generate_log_decay_coefs <- function() {
    a <- rnorm(1, mean = coef_estimates[[1]], sd = coef_se[[1]])
    b <- rnorm(1, mean = coef_estimates[[2]], sd = coef_se[[2]])
    return(c(a, b))
  }
  
  log_linear_decay <- function(t, a, b) {
    log_odn <- a + b * t
    pmax(exp(log_odn), 0.001)
  }
  
  sigmoid <- function(t, midpoint, slope = 2) {
    1 / (1 + exp(-slope * (t - midpoint)))
  }
  
  logistic <- function(t, midpoint, slope = 6) {
    1 / (1 + exp(-slope * (t - midpoint)))
  }
  
  exp_rebound <- function(start_value, t, rebound_rate, fail_time) {
    rebound <- start_value * exp(rebound_rate * (t - fail_time))
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
      times <- seq(0, fail_time + 3 * time_interval, by = time_interval)
      censor_time <- max(times)
    } else if (dropout_flags[i]) {
      times <- seq(0, dropout_times[i], by = time_interval)
      censor_time <- dropout_times[i]
    } else {
      times <- full_time_points
      censor_time <- max_follow_up
    }
    
    follow_up_years[i] <- max(times)
    coefs <- generate_log_decay_coefs()
    a <- coefs[1]
    b <- coefs[2]
    decay_params[i, ] <- coefs
    
    if (has_failed) {
      fail_idx <- which(times >= fail_time)[1]
      treatment_failure_vector[i] <- TRUE
      rebound_rate <- runif(1, min = 0.01, max = 0.05)
      rebound_rate_vector[i] <- rebound_rate
      
      decay_vals <- log_linear_decay(times, a, b)
      rebound_start <- decay_vals[fail_idx]
      rebound_vals <- exp_rebound(rebound_start, times, rebound_rate, fail_time)
      
      transition_weights <- if (rebound_model == "sigmoid") {
        sigmoid(times, fail_time)
      } else {
        logistic(times, fail_time)
      }
      
      combined_odn <- (1 - transition_weights) * decay_vals + transition_weights * rebound_vals
    } else {
      combined_odn <- log_linear_decay(times, a, b)
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

result <- simulate_ODn_decay2(
  coef_estimates = coef_estimates,
  coef_se = coef_se,
  n_individuals = 10000,
  baseline_mean = 3.47,
  baseline_sd = 1.55,
  sigma_0 = -0.01469,
  slope_sigma = 0.14513,
  baseline_noise = 0.1,
  fraction = 0.15,
  max_follow_up = 10,
  time_interval = 0.5,
  noise_model = 1, # 1 = heteroskedastic, 2 = realistic noise
  failure_prob = 0.25,
  dropout_prob = 0.1,
  rebound_model = "sigmoid"  # or "logistic"
)
decay_data <- result$decay_data
model_parameters <- result$model_parameters

#' Select the Best-Matching Log-Linear Decay Model for a Test Individual
#'
#' Given a test individual's longitudinal ODn data and a parameter set of simulated models
#' (excluding those with treatment failure), this function finds the best-fitting decay model
#' using root mean squared error (RMSE), mean absolute error (MAE), and a simple model complexity metric.
#'
#' @param test_data A data frame containing the observed data for a single individual.
#'   It must have a column named `record_id`, followed by columns representing time and ODn values
#'   (usually named `time` and `value` respectively).
#' @param param_sim_data A data frame containing simulated model parameters for multiple individuals.
#'   Must include columns: `id`, `baseline`, `a`, `b`, and `treatment_failure`.
#'
#' @return A data frame with the best-matching model parameters, including:
#' \describe{
#'   \item{id}{ID of the selected simulated individual.}
#'   \item{baseline}{Baseline ODn value.}
#'   \item{a, b}{Log-linear decay coefficients.}
#'   \item{rmse}{Root mean squared error between model predictions and observed ODn.}
#'   \item{mae}{Mean absolute error.}
#'   \item{record_id}{ID of the test individual.}
#' }
#'
#' @details
#' The function matches the observed ODn trajectory of a test individual against a library
#' of simulated trajectories defined by log-linear decay parameters (`a`, `b`), using the following criteria:
#' \enumerate{
#'   \item Lowest RMSE
#'   \item Lowest MAE (used to break ties)
#'   \item Lowest model complexity, defined as \code{abs(a) + abs(b)} (further tie-breaker)
#'   \item Random choice if still tied
#' }
#'
#' Simulated models with treatment failure are excluded. Infinite or NaN errors in predictions
#' are filtered out prior to selection.
#'
#' @importFrom dplyr %>%, filter, select
#'
#' @examples
#' # Example usage (assuming simulated parameter data and test data are available)
#' best_model <- best_model_choice2(test_data = test_individual_df,
#'                                  param_sim_data = simulation_df)
#' print(best_model)
#'
#' @export
best_model_choice2 <- function(test_data, param_sim_data) {
  library(dplyr)
  
  # Get the record_id from the first row of test_data (assumes single individual)
  record_id <- test_data$record_id[1]
  
  # Filter out individuals who had treatment failure (i.e., viral rebound)
  param_sim_data <- param_sim_data %>% filter(!treatment_failure)
  
  if (nrow(param_sim_data) == 0) {
    stop("No eligible models available: all individuals had viral rebound.")
  }
  
  # Extract time and ODn values from test data
  time_vec <- test_data[[2]]
  ODn_vec <- test_data[[3]]
  n_timepoints <- length(time_vec)
  n_models <- nrow(param_sim_data)
  
  # Create prediction matrix using log-linear model
  time_matrix <- matrix(rep(time_vec, each = n_models), nrow = n_models)
  a_matrix <- matrix(param_sim_data$a, nrow = n_models, ncol = n_timepoints)
  b_matrix <- matrix(param_sim_data$b, nrow = n_models, ncol = n_timepoints)
  
  pred_matrix <- exp(a_matrix + b_matrix * time_matrix)
  pred_matrix <- pmax(pred_matrix, 0.001)
  
  # Actual observed ODn values replicated
  obs_matrix <- matrix(rep(ODn_vec, each = n_models), nrow = n_models)
  
  # Error matrices
  sq_error_matrix <- (obs_matrix - pred_matrix)^2
  abs_error_matrix <- abs(obs_matrix - pred_matrix)
  
  # Filter out rows with infinite errors
  finite_rows <- is.finite(rowSums(sq_error_matrix))
  param_sim_data <- param_sim_data[finite_rows, , drop = FALSE]
  sq_error_matrix <- sq_error_matrix[finite_rows, , drop = FALSE]
  abs_error_matrix <- abs_error_matrix[finite_rows, , drop = FALSE]
  
  if (nrow(param_sim_data) == 0) {
    stop("All filtered models resulted in invalid predictions.")
  }
  
  # Performance metrics
  param_sim_data$rmse <- sqrt(rowMeans(sq_error_matrix))
  param_sim_data$mae <- rowMeans(abs_error_matrix)
  param_sim_data$complexity <- abs(param_sim_data$a) + abs(param_sim_data$b)
  
  # Select best model
  dt_min <- param_sim_data %>% filter(rmse == min(rmse))
  if (nrow(dt_min) > 1) {
    dt_min <- dt_min %>% filter(mae == min(mae))
  }
  if (nrow(dt_min) > 1) {
    dt_min <- dt_min %>% filter(complexity == min(complexity))
  }
  if (nrow(dt_min) > 1) {
    set.seed(11)
    dt_min <- dt_min[sample(1:nrow(dt_min), 1), ]
  }
  
  # Add record_id to output
  dt_min$record_id <- record_id
  
  return(dt_min %>% select(id, baseline, a, b, rmse, mae, record_id))
}

# library(dplyr)
library(purrr)

# Filter once outside the loop
test_data_filtered <- test_data %>% 
  filter(is.na(peak_visit)) %>% 
  group_split(record_id)

# Apply best_model_choice efficiently using map
results_list <- map(test_data_filtered, function(individual_data) {
  best_model_choice2(
    test_data = individual_data,
    param_sim_data = model_parameters
  )
})

# Combine all individual results into a single dataframe
results <- results_df <- bind_rows(results_list)

#' Compare an Observed ODn Value to a Log-Linear Model Prediction
#'
#' This function compares an observed ODn value at a given time point to a model-predicted value
#' using a log-linear decay model. It calculates the z-statistic and corresponding one-sided
#' p-value for the hypothesis that the observed ODn is significantly greater than expected.
#'
#' @param a Numeric. The intercept parameter from the log-linear decay model.
#' @param b Numeric. The slope parameter from the log-linear decay model.
#' @param baseline Numeric. The baseline ODn value for the individual.
#' @param t Numeric. The time point (in years) at which the observed value is compared.
#' @param y_ODn Numeric. The observed ODn value at time `t`.
#' @param sigma_ODn Numeric. The standard deviation of the ODn prediction at time `t`.
#' @param sigma_y_ODn Numeric. The measurement error (standard deviation) of the observed ODn.
#' @param id Identifier for the individual (can be numeric or character).
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{id}{Identifier for the individual.}
#'   \item{value}{The observed ODn value.}
#'   \item{predicted}{The expected ODn value under the model.}
#'   \item{z_stat}{The z-statistic comparing observed vs. predicted values.}
#'   \item{p_value}{The one-sided p-value for \eqn{H_1: y_{ODn} > \hat{y}}.}
#'   \item{flagged}{Logical indicator: \code{TRUE} if the result is statistically significant (p < 0.05).}
#' }
#'
#' @details
#' The predicted ODn is computed using the formula:
#' \deqn{\hat{y} = \max(baseline \cdot \exp(-(a + b \cdot t) \cdot t), 0.001)}
#' 
#' The pooled standard deviation assumes independent error sources and is computed as:
#' \deqn{\sigma_{pooled} = \sqrt{\sigma_{ODn}^2 + \sigma_{y}^2}}
#'
#' The z-statistic and one-sided p-value test whether the observed value is significantly
#' greater than expected under the null model. This is useful for identifying individuals
#' who may have ODn trajectories inconsistent with typical decay patterns.
#'
#' @importFrom tibble tibble
#' @importFrom stats pnorm
#'
#' @examples
#' compare_value_with_others2(
#'   a = 0.3, b = -0.1, baseline = 2.0,
#'   t = 2, y_ODn = 1.5,
#'   sigma_ODn = 0.2, sigma_y_ODn = 0.15,
#'   id = "patient_001"
#' )
#'
#' @export

compare_value_with_others2 <- function(a, b, baseline, t, y_ODn, sigma_ODn, sigma_y_ODn, id) {
  # Compute predicted ODn value from log-linear model
  y_hat <- pmax(baseline * exp(-(a + b * t) * t), 0.001)
  
  # Compute pooled standard deviation
  sigma_pooled <- sqrt(sigma_ODn^2 + sigma_y_ODn^2)
  
  # Compute z-statistic
  z_stat <- (y_ODn - y_hat) / sigma_pooled
  
  # One-sided p-value for H1: y_ODn > y_hat
  p_value <- 1 - pnorm(z_stat)
  
  # Significance flag
  flagged <- p_value < 0.05
  
  # Return results as a tibble
  tibble::tibble(
    id = id,
    value = y_ODn,
    predicted = y_hat,
    z_stat = z_stat,
    p_value = p_value,
    flagged = flagged
  )
}

# Step 1: Filter last visits (post-rebound)
test_data_last_visit <- test_data %>% 
  filter(!is.na(peak_visit)) %>% 
  group_by(record_id) %>% 
  slice_tail(n = 1) %>% 
  ungroup()

# Step 2: Join with model parameters
data_combined <- test_data_last_visit %>% 
  left_join(results, by = c("record_id"))

# Step 3: Efficiently compute pooled SDs and p-values
results_final <- pmap_dfr(
  list(
    split(data_combined, seq_len(nrow(data_combined))),
    data_combined$time_vec,
    data_combined$ODn_vec,
    map_dbl(data_combined$record_id, ~ sd(test_data$ODn_vec[test_data$record_id == .x & is.na(test_data$peak_visit)])),
    noise$coefficients[1, 1] + noise$coefficients[2, 1] * data_combined$time_vec,
    data_combined$record_id
  ),
  ~ compare_value_with_others(..1, ..2, ..3, ..4, ..5, ..6)
)

dt05 <- as.data.frame(results_final) %>%
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
  dplyr::select(strata, flagged) %>%
  tbl_summary(by = strata
  )
dt06

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