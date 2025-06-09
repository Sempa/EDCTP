library(readxl)
library(tidyverse)
library(gtsummary)
library(dplyr)
library(ggplot2)
library(minpack.lm)
library(nlme)
library(splines)
library(cowplot)
library(distributions3)
library("gghighlight")
library(brolgar)
library(glmmTMB)
library(purrr)

sediaData_full <- read_csv("data/JHU/CEPHIA - JHU LAg-Avidity Data.csv")

sediaData <- sediaData_full %>%
  filter(!is.na(days_since_eddi)) %>%
  filter(!is.na(viral_load)) %>%
  filter(on_treatment == 'TRUE') %>%
  mutate(subject_eddi = paste(subject_label, days_since_eddi, sep = ' ')) %>%
  distinct(subject_eddi, .keep_all = T) %>%
  group_by(subject_label) %>%
  mutate(n_visits = 1:length(subject_label),
         sup_visits = ifelse(viral_load<1000,1,2)) %>%
  mutate(tot_visits = max(n_visits),
         tot_supp_visits = sum(sup_visits)) %>%
  mutate(fully_supp = ifelse(tot_visits == tot_supp_visits, 1, 0)) %>%
  ungroup() %>%
  filter(tot_visits>2) %>%
  # filter(subject_label %in% x) %>%
  dplyr::select(subject_label_blinded = subject_label, sex, viral_load, days_since_eddi, ODn, n_visits) %>%
  arrange(subject_label_blinded, days_since_eddi)
dt <- sediaData %>%
  mutate(x= ifelse(viral_load > 999, 1,0)) %>%
  group_by(subject_label_blinded) %>%
  mutate(flag = max(x, na.rm = T)) %>%
  filter(flag == 0) %>%
  arrange(subject_label_blinded, days_since_eddi)

# distinct(subject_label, .keep_all = T)
baseline_ODn_data <- sediaData %>%
  mutate(flag = ifelse(viral_load > 10000 & n_visits==1, 1,0)) %>%
  filter(flag==1)
cephia_samples <- read_csv("Sempa_final_pull_with_results.csv") %>%
  mutate(vl = ifelse(`Viral Load at Draw` == "<40", "40", ifelse(`Viral Load at Draw` == "NULL", "", `Viral Load at Draw`))) %>%
  mutate(
    logvl = log(as.numeric(vl), 10),
    time_on_trt = -1 * `days since tx start`,
    sedia_ODn = `Sedia LAg Odn screen`,
    RaceEthnicity = `Race/Ethnicity`,
    viral_load = as.numeric(vl),
    test_date = visit_date,
    days_since_eddi = `Days from EDDI to draw`
  ) %>%
  arrange(subject_label_blinded, time_on_trt) %>%
  group_by(subject_label_blinded) %>%
  mutate(
    visits = 1:length(subject_label_blinded),
    unsupressed_visits = ifelse(viral_load > 999, 1, 0)
  ) %>%
  mutate(
    visits = ifelse(unsupressed_visits == 1, visits, NA),
    baseline_visit = ifelse(unsupressed_visits == 1, max(visits, na.rm = T), 0)
  ) %>%
  mutate(baseline_visit = ifelse((subject_label_blinded == 35329295 | subject_label_blinded == 52033420), 1,
                                 ifelse(subject_label_blinded == 54382839, 5,
                                        ifelse(subject_label_blinded == 64577680, 13, baseline_visit)
                                 )
  )) %>%
  mutate(
    baseline_visit = max(baseline_visit),
    visits = 1:length(subject_label_blinded)
  ) %>%
  filter(visits >= baseline_visit) %>%
  mutate(
    `days since tx start` = `days since tx start` * -1,
    subject_label_blinded = as.character(subject_label_blinded),
    to_peak = NA
  ) %>%
  select(subject_label_blinded, days_since_eddi, `days since tx start`,
         Sex = BiologicalSex,# Subtype,
         Age = `Age at Draw`, test_date, sedia_ODn, viral_load, visits, Group, to_peak
  )
full_dataset <- bind_rows(
  cephia_samples %>%
    mutate(cohort = "cephia"),
  read_delim("data/full_africos_data_with_ODn.csv",
             delim = ";", escape_double = FALSE, trim_ws = TRUE
  ) %>% # africos_pts <-
    filter(exclude %in% c(1, 2)) %>%
    mutate(
      days_since_eddi = NA,
      `days since tx start` = NA,
      Group = ifelse(exclude == 1, "early suppressions",
                     ifelse(exclude == 2, "suppression failures", NA)
      ),
      Group = ifelse(is.na(Group), "suppression failures", Group),
      Age = ifelse(Age == ".", "31.6", Age)
    ) %>%
    select(
      subject_label_blinded = subjid, days_since_eddi, `days since tx start`, Sex, Age,
      test_date = Date_Collected, sedia_ODn = ODn, viral_load = vl,
      visits = Visit, Group, exclude, to_peak
    ) %>%
    left_join(
      readxl::read_excel("data/AFRICOS_Sempa Data Pull_24Jul23.xlsx",
                         sheet = "Sheet1"
      ) %>%
        mutate(id = paste(Date_Collected, subjid, sep = "_")) %>%
        group_by(subjid) %>%
        mutate(
          switch_date = as.Date(as.numeric(ifelse(`ARV Start Date` == ".", "", `ARV Start Date`)),
                                origin = "1900-01-01"
          ),
          art_start = as.Date(min(as.numeric(ifelse(`ARV Start Date` == ".", "", `ARV Start Date`)), na.rm = T),
                              origin = "1900-01-01"
          )
        ) %>%
        ungroup() %>%
        dplyr::select(subject_label_blinded = subjid, art_start) %>%
        distinct(subject_label_blinded, .keep_all = T),
      by = "subject_label_blinded"
    ) %>%
    mutate(
      `days since tx start` = as.numeric(as.Date(test_date, origin = "1900-01-01") - art_start),
      cohort = "africos",
      # subject_label_blinded = as.double(cur_group_id()),
      Age = as.numeric(Age),
      viral_load = as.numeric(viral_load)
    ) %>%
    select(
      subject_label_blinded, days_since_eddi,
      `days since tx start`, Sex, Age,
      test_date, sedia_ODn, viral_load, visits, Group, cohort, to_peak
    )
)  %>%
  mutate(years_since_tx_start = `days since tx start`/365.25)

sedia_eddi_data <- bind_rows(
  dt %>%
    dplyr::select(subject_label_blinded, sex, days_since_eddi, sedia_ODn = ODn),
  cephia_samples %>%
    filter(Group == 'early suppressions') %>%
    dplyr::select(subject_label_blinded, sex = Sex, days_since_eddi, sedia_ODn)
)

baseline_ODn_table <- full_dataset %>%
  group_by(subject_label_blinded) %>%
  mutate(n_visit = 1:length(subject_label_blinded)) %>%
  ungroup() %>%
  filter(n_visit == 1)
mean(c(baseline_ODn_table$sedia_ODn, baseline_ODn_data$ODn)) ## baseline mean
sd(c(baseline_ODn_table$sedia_ODn, baseline_ODn_data$ODn))

# model_eddi <- nlme::lme(sedia_ODn ~ poly(days_since_eddi, 2), 
#                   data = sedia_eddi_data, 
#                   random = ~ days_since_eddi|subject_label_blinded,
#                   na.action = na.omit)
# summary(model_eddi)
# ploy_model <- nlme::lme(sedia_ODn ~ poly(years_since_tx_start, 2), 
#                         data = full_dataset %>%
#                           filter(Group == 'early suppressions'), 
#                         random = ~years_since_tx_start|subject_label_blinded,
#                         na.action = na.omit)
# Use the coefficients as starting values
poly_model <- glmmTMB(
  sedia_ODn ~ years_since_tx_start + (1 | subject_label_blinded),
  family = gaussian(link = "log"),
  data = full_dataset,
  start = list(beta = c(0, 1))
)
coefs <- summary(poly_model)$coefficients$cond
coef_estimates <- coefs[, "Estimate"]
coef_se <- coefs[, "Std. Error"]
summary(poly_model)

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
###########################################################
##getting data to use for testing - CEPHIA EDCTP data
###########################################################
dt02 <- read_csv("data/20180410-EP-LAgSedia-Generic.csv") %>%
  filter(on_treatment == T) %>%
  dplyr::select(subject_label_blinded, visit_date, test_date, art_initiation_date, 
                days_since_eddi, viral_load, sedia_ODn = result...72) %>% 
  arrange(subject_label_blinded, test_date) %>%
  # filter(as.character(art_initiation_date) != '') %>%
  # filter(visit_date >= art_initiation_date) %>%
  mutate(years_since_ART_start = as.numeric(visit_date-art_initiation_date)/365.25) %>%
  arrange(subject_label_blinded, visit_date) %>%
  group_by(subject_label_blinded) %>%
  mutate(x = 1:length(subject_label_blinded)) %>%
  mutate(flag = max(x)) %>%
  ungroup() %>%
  filter(flag >2) %>%
  group_by(subject_label_blinded) %>%
  mutate(x = 1:length(subject_label_blinded),
         flag = max(x)) %>%
  mutate(peak = ifelse(x==flag, 1, NA)) %>%
  ungroup() %>%
  mutate(strata = 'early suppressions - cephia') %>%
  dplyr::select(subject_label_blinded, strata, years_since_ART_start,sedia_ODn, peak)

dt04 <- full_dataset %>%
  filter(Group == 'suppression failures') %>%
  group_by(subject_label_blinded) %>%
  mutate(x = 1:length(subject_label_blinded),
         flag = max(x)) %>%
  mutate(peak = ifelse(x==flag, 1, NA)) %>%
  ungroup() %>%
  filter(flag >2) %>%
  mutate(years_since_ART_start = `days since tx start`/365.25) %>%
  dplyr::select(subject_label_blinded, strata = Group, years_since_ART_start, sedia_ODn, peak)
test_data <- bind_rows(
  dt02 %>%
    mutate(subject_label_blinded = as.character(subject_label_blinded)), dt04
) %>%
  mutate(id = as.numeric(as.factor(subject_label_blinded))) %>%
  dplyr::select(record_id = id, time_vec = years_since_ART_start, 
                ODn_vec = sedia_ODn, peak_visit = peak)

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
results <- results_df %>%
  mutate(record_id = ids)

###########################################################
###sigma ODn
############################################################
controls_blinded_sedia <-  read_csv("data/20180410-EP-LAgSedia-Generic.csv") %>%
  mutate(sedia_ODn = `result...15`) %>%
  filter(visit_id %in% c(21773, 21783, 21785)) %>%
  mutate(specimen = case_when(
    visit_id == 21773 ~ "BC-1",
    visit_id == 21783 ~ "BC-2",
    visit_id == 21785 ~ "BC-3"
  )) %>%
  select(specimen, test_date,
         specimen_blinded_label = specimen_label,
         sedia_final_ODn = `result...15`,
         sedia_final_ODn_method = method
  )
sedia_distribution_blinded_controls <- controls_blinded_sedia %>%
  group_by(specimen) %>%
  summarise(
    `mean Sedia ODn` = round(mean(sedia_final_ODn), 3),
    `sigma Sedia ODn` = round(sd(sedia_final_ODn), 3)
  ) %>%
  mutate(`sigma over mean ODn` = `sigma Sedia ODn` / `mean Sedia ODn`)

sd_sedia_ODn_vs_ODn_bc <- sedia_distribution_blinded_controls %>%
  ggplot(aes(x = `mean Sedia ODn`, y = `sigma Sedia ODn`)) +
  geom_point() +
  expand_limits(x = 0, y = 0) +
  geom_smooth(method = lm, size = 1.5, se = FALSE)

sd_sedia_ODn_vs_ODn_bc

noise <- summary(glm(sd_Sedia_ODn ~ mean_Sedia_ODn, data = sedia_distribution_blinded_controls %>%
                       mutate(sd_Sedia_ODn = `sigma Sedia ODn`, mean_Sedia_ODn = `mean Sedia ODn`)))
######################################################################
#How about a new time point?
######################################################################
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
