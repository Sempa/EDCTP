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

# 1. Filter out individuals with insufficient data (e.g., <2 points)
filtered_subjects <- full_dataset %>%
  group_by(subject_label_blinded) %>%
  filter(n() >= 2) %>%
  # filter(Group == 'early suppressions') %>%
  ungroup()

# 2. Fit glm with log-link on each individual
individual_models <- filtered_subjects %>%
  group_split(subject_label_blinded) %>%
  map_df(function(df) {
    # Fit glm on log scale
    fit <- try(glm(log(sedia_ODn) ~ years_since_tx_start, data = df), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
      return(NULL)
    }
    
    coefs <- coef(fit)
    tibble(
      subject_label_blinded = df$subject_label_blinded[1],
      intercept = coefs[1],
      slope = coefs[2]
    )
  })

# 3. Calculate mean and SD of intercept and slope
intercept_mean <- mean(individual_models$intercept, na.rm = TRUE)
intercept_sd   <- sd(individual_models$intercept, na.rm = TRUE)
slope_mean     <- mean(individual_models$slope, na.rm = TRUE)
slope_sd       <- sd(individual_models$slope, na.rm = TRUE)
coef_estimates <- list(intercept_mean, slope_mean)
coef_se <- list(intercept_sd, slope_sd)

generate_model_parameters <- function(
    coef_estimates, coef_se, 
    n_individuals, failure_rate = 0.1) {
  # browser()
  set.seed(123)
  
  library(truncnorm)
  
  a0_values <- truncnorm::rtruncnorm(n_individuals, a = 0.03, b = 7.4,
                                     mean = coef_estimates[[1]], sd = coef_se[[1]])
  lambda_i_values <- rnorm(n_individuals, mean = coef_estimates[[2]], sd = coef_se[[2]])
  fail_times <- log(runif(n_individuals, min = 0, max = 1)) / (-1 * failure_rate)
  gamma_values <- runif(n_individuals, min = 0.01, max = 0.05)
  
  model_parameters <- tibble::tibble(
    id = 1:n_individuals,
    a = a0_values,
    b = lambda_i_values,
    lambda_i = lambda_i_values,
    gamma = gamma_values,
    failure_time = fail_times
  )
  
  return(model_parameters)
}

model_parameters <- generate_model_parameters(
  coef_estimates, 
  coef_se, 
  n_individuals = 1000)

summary(model_parameters$failure_time)

generate_decay_data <- function(model_parameters, 
                                max_follow_up = 10,
                                time_interval = 0.5,
                                noise_model = 1,
                                sigma_0 = 0.1,
                                slope_sigma = 0.05,
                                baseline_noise = 0.1,
                                fraction = 0.1) {
  
  full_time_points <- seq(0, max_follow_up, by = time_interval)
  compute_noise_sd1 <- function(odn) pmax(sigma_0 + slope_sigma * odn, 0.01)
  compute_noise_sd2 <- function(odn) pmax(baseline_noise + fraction * odn, 0.01)
  
  decay_data_list <- vector("list", nrow(model_parameters))
  
  for (i in seq_len(nrow(model_parameters))) {
    a0 <- model_parameters$a[i]
    lambda_i <- model_parameters$lambda_i[i]
    gamma <- model_parameters$gamma[i]
    failure_time <- model_parameters$failure_time[i]
    times <- full_time_points
    
    odn_true <- ifelse(times < failure_time,
                       a0 * exp(-lambda_i * times),
                       pmin(a0 * (1 - (1 - exp(-lambda_i * failure_time)) * 
                                    exp(-gamma * (times))), a0))
    
    noise_sd <- if (noise_model == 1) compute_noise_sd1(odn_true) else compute_noise_sd2(odn_true)
    odn_noisy <- odn_true
    if (length(odn_noisy) > 1) {
      odn_noisy[-1] <- pmin(
        pmax(rnorm(length(times) - 1, mean = odn_true[-1], sd = noise_sd[-1]), 0.001),
        a0
      )
    }
    
    decay_data_list[[i]] <- tibble::tibble(
      individual = model_parameters$id[i],
      time = times,
      value = odn_noisy,
      failure_time = failure_time
    )
  }
  
  decay_data <- dplyr::bind_rows(decay_data_list)
  return(decay_data)
}

# Generate decay data using those parameters:
decay_data <- generate_decay_data(
  model_parameters,
  sigma_0 = 0.1,
  slope_sigma = 0.05,
  baseline_noise = 0.2,
  fraction = 0.3,
  max_follow_up = 10,
  time_interval = 0.5,
  noise_model = 1
)

# Call simulate_ODn_3
# results <- simulate_ODn_3(
#   coef_estimates = coef_estimates,
#   coef_se = coef_se,
#   n_individuals = 1000,
#   sigma_0 = 0.2,
#   slope_sigma = 0.1,
#   baseline_noise = 0.1,
#   fraction = 0.1,
#   max_follow_up = 10,
#   time_interval = 0.5,
#   noise_model = 1
# )

# Access outputs
decay_data <- results$decay_data
model_parameters <- results$model_parameters

n_individuals <- 1000
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

ggplot_plots

detect_ODn_upticks <- function(
    decay_data, 
    sd_option = c("fixed", "rolling_window", "rolling_all"), 
    z_threshold = 1.96) {
  
  library(dplyr)
  library(tidyr)
  
  sd_option <- match.arg(sd_option)
  fixed_sd <- 1.85  # Sedia control-derived fixed SD
  
  decay_data <- decay_data %>%
    arrange(individual, time) %>%
    group_by(individual) %>%
    mutate(
      z_score = sapply(seq_along(value), function(i) {
        if (i < 3) return(NA_real_)
        prev_vals <- value[max(1, i - 4):(i - 1)]
        if (length(prev_vals) < 1) return(NA_real_)
        
        mean_prev <- mean(prev_vals)
        
        sd_val <- switch(sd_option,
                         fixed = fixed_sd,
                         rolling_window = if (length(prev_vals) > 1) sd(prev_vals) else NA_real_,
                         rolling_all = if (i > 2) sd(value[1:(i - 1)]) else NA_real_)
        
        if (is.na(sd_val) || sd_val == 0) return(NA_real_)
        
        (value[i] - mean_prev) / sd_val
      }),
      ODn_uptick_flag = z_score > z_threshold
    ) %>%
    mutate(
      uptick_time = ifelse(row_number() == which(ODn_uptick_flag)[1], time, NA_real_)
    ) %>%
    fill(uptick_time, .direction = "down") %>%
    mutate(
      uptick_time = ifelse(row_number() == 1, NA_real_, uptick_time),
      time_since_failure = ifelse(ODn_uptick_flag,
                                  ifelse(!is.na(failure_time), time - failure_time, Inf),
                                  NA_real_)
    ) %>%
    ungroup()
  
  return(decay_data)
}

# updated_decay_data <- detect_ODn_upticks(decay_data, model_parameters, sd_option = "rolling_window", z_threshold = 1.96)


# Example usage
# Assuming `decay_data` is already defined and loaded:
# Detect upticks using fixed SD and threshold of 1.96
updated_decay_data <- detect_ODn_upticks(decay_data, sd_option = "fixed", z_threshold = 1.96)
table(updated_decay_data$ODn_uptick_flag)
# Use rolling window SD
updated_decay_data2 <- detect_ODn_upticks(decay_data, sd_option = "rolling_window", z_threshold = 1.96)
table(updated_decay_data2$ODn_uptick_flag)
# Use cumulative rolling SD
updated_decay_data3 <- detect_ODn_upticks(decay_data, sd_option = "rolling_all", z_threshold = 1.96)
table(updated_decay_data3$ODn_uptick_flag)

# Summarize at individual level
summary_df <- updated_decay_data2 %>%
  group_by(individual) %>%
  summarise(
    has_uptick = any(ODn_uptick_flag, na.rm = TRUE),
    uptick_time = if (any(ODn_uptick_flag, na.rm = TRUE)) {
      min(time[ODn_uptick_flag], na.rm = TRUE)
    } else {
      NA_real_
    },
    failure_time = unique(failure_time),
    time_diff = if (!is.na(uptick_time)) uptick_time - failure_time else NA_real_,
    .groups = "drop"
  )

# Counts
n_without_uptick <- sum(!summary_df$has_uptick)
n_with_uptick <- sum(summary_df$has_uptick) # & (summary_df$failure_time >1 & summary_df$failure_time <= 9)
n_fail_lt_10_with_uptick <- sum(summary_df$has_uptick & summary_df$failure_time < 10)
n_fail_gt_10_with_uptick <- sum(summary_df$has_uptick & summary_df$failure_time >= 10)
delays_to_detection_of_rebound <- quantile(summary_df$time_diff, na.rm = TRUE)

# View or print results
print(paste("Number without ODn uptick:", n_without_uptick))
print(paste("Number with ODn uptick:", n_with_uptick))
print(paste("Number with ODn uptick and failure_time < 10:", n_fail_lt_10_with_uptick))
print(paste("Number with ODn uptick and failure_time >= 10:", n_fail_gt_10_with_uptick))
print(paste("Distribution of delays to detection (uptick time-failure time)",delays_to_detection_of_rebound))

# Add these summary metrics to the individual-level data
summary_df <- summary_df %>%
  mutate(
    uptick_category = case_when(
      has_uptick & failure_time < 10 ~ "Uptick with failure < 10",
      has_uptick & failure_time >= 10 ~ "Uptick with failure ≥ 10",
      !has_uptick ~ "No uptick"
    )
  )

  