poly_model <- glmmTMB(
  sedia_ODn ~ years_since_tx_start + (1 | subject_label_blinded),
  family = gaussian(link = "log"),
  data = full_dataset,
  start = list(beta = c(0, 1))
)
coefs <- summary(poly_model)$coefficients$cond
coef_estimates <- coefs[, "Estimate"]
coef_se <- coefs[, "Std. Error"]

set.seed(123)

# Number of individuals
n_individuals <- 10000

# Time points (3-month intervals up to 10 years)
all_time_points <- seq(0, 10, by = 0.25)

# Baseline ODn values
baseline_mean <- 3.47
baseline_sd <- 1.55
baselines <- rtruncnorm(n_individuals, mean = baseline_mean, sd = baseline_sd, a = 0.03, b = 7.4)

# Polynomial coefficients and standard errors
coef_estimates <- list(a = 1.818730, b = -8.97018)
coef_se <- list(a = 0.2, b = 0.5)

# Generate decay rates from polynomial
generate_decay_rate <- function(t) {
  a <- rnorm(1, mean = coef_estimates$a, sd = coef_se$a)
  b <- rnorm(1, mean = coef_estimates$b, sd = coef_se$b)
  decay_rate <- pmax(a + b * t, 0.05)
  return(list(decay_rate = decay_rate, a = a, b = b))
  library(dplyr) 
  library(truncnorm)
  
  set.seed(123)
  
  # Number of individuals
  n_individuals <- 10000
  
  # Time points (3-month intervals up to 10 years)
  all_time_points <- seq(0, 10, by = 0.25)
  
  # Baseline ODn values
  baseline_mean <- 3.47
  baseline_sd <- 1.55
  baselines <- rtruncnorm(n_individuals, mean = baseline_mean, sd = baseline_sd, a = 0.03, b = 7.4)
  
  # Polynomial coefficients and standard errors
  coef_estimates <- list(a = 1.818730, b = -8.97018)
  coef_se <- list(a = 0.2, b = 0.5)
  
  # Generate decay rates from polynomial
  generate_decay_rate <- function(t) {
    a <- rnorm(1, mean = coef_estimates$a, sd = coef_se$a)
    b <- rnorm(1, mean = coef_estimates$b, sd = coef_se$b)
    decay_rate <- pmax(a + b * t, 0.05)
    return(list(decay_rate = decay_rate, a = a, b = b))
  }
  
  # Exponential decay function
  exp_decay <- function(t, baseline, decay_rate) {
    pmax(baseline * exp(-decay_rate * t), 0.001)
  }
  
  # Noise model 1: Regression-based
  compute_noise_sd1 <- function(odn) {
    pmax(-0.01469 + 0.14513 * odn, 0.01)
  }
  
  # Noise model 2: Base + proportional
  compute_noise_sd2 <- function(odn, base_sd = 0.1, prop_sd = 0.05) {
    pmax(base_sd + prop_sd * odn, 0.01)
  }
  
  # Choose noise model (1 or 2)
  noise_model <- 2
  
  # Skewed follow-up times (mostly 10 years)
  follow_up_years <- round(qbeta(runif(n_individuals), 5, 1) * 10, 1)
  follow_up_years[follow_up_years < 1] <- 1
  
  # Simulate treatment failure
  fail_rate <- 4.5 / 100
  expected_failures <- rbinom(n_individuals, 1, prob = 1 - exp(-fail_rate * follow_up_years))
  failure_time <- rexp(n_individuals, rate = fail_rate)
  failure_time <- pmin(failure_time, follow_up_years)
  failure_time[expected_failures == 0] <- NA
  
  # Simulate longitudinal data
  decay_data <- do.call(rbind, lapply(1:n_individuals, function(i) {
    t_i <- all_time_points[all_time_points <= follow_up_years[i]]
    baseline_i <- baselines[i]
    
    failed_i <- !is.na(failure_time[i])
    fail_time_i <- failure_time[i]
    
    if (!failed_i) {
      decay_i <- generate_decay_rate(t_i)
      expected_odn <- exp_decay(t_i, baseline_i, decay_i$decay_rate)
    } else {
      t_pre <- t_i[t_i <= fail_time_i]
      t_post <- t_i[t_i > fail_time_i]
      
      decay_pre <- generate_decay_rate(t_pre)
      expected_pre <- exp_decay(t_pre, baseline_i, decay_pre$decay_rate)
      
      last_pre_index <- length(t_pre)
      t_last <- t_pre[last_pre_index]
      
      decay_post <- generate_decay_rate(t_post - fail_time_i)
      baseline_post <- expected_pre[last_pre_index]
      expected_post <- exp_decay(t_post - fail_time_i, baseline_post, decay_post$decay_rate)
      
      # Replace last pre-failure ODn with failure trajectory
      expected_pre[last_pre_index] <- exp_decay(0, baseline_post, decay_post$decay_rate)[1]
      
      expected_odn <- c(expected_pre, expected_post)
      t_i <- c(t_pre, t_post)
    }
    
    noise_sd <- switch(as.character(noise_model),
                       `1` = compute_noise_sd1(expected_odn),
                       `2` = compute_noise_sd2(expected_odn))
    
    observed_odn <- pmax(rnorm(length(t_i), mean = expected_odn, sd = noise_sd), 0.001)
    
    data.frame(
      individual = i,
      time = t_i,
      expected_odn = expected_odn,
      observed_odn = observed_odn,
      noise_sd = noise_sd,
      failed = failed_i,
      failure_time = ifelse(failed_i, fail_time_i, NA)
    )
  }))
  
  # Store model parameters
  decay_params <- t(sapply(1:n_individuals, function(i) {
    t_i <- all_time_points[all_time_points <= follow_up_years[i]]
    decay_i <- generate_decay_rate(t_i)
    c(a = decay_i$a, b = decay_i$b)
  }))
  
  model_parameters <- data.frame(
    id = 1:n_individuals,
    baseline = baselines,
    follow_up_years = follow_up_years,
    failure_time = failure_time,
    failed = !is.na(failure_time),
    a = decay_params[, "a"],
    b = decay_params[, "b"]
  )
  
