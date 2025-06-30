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
  filter(Group == 'early suppressions') %>%
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

simulate_ODn_3 <- function(
    coef_estimates, coef_se,
    n_individuals,
    sigma_0, slope_sigma,
    baseline_noise, fraction,
    failure_prob = 0.2,
    max_follow_up = 10,
    time_interval = 0.5,
    noise_model = 1
) {
  # library(dplyr)
  
  # Generate time vector
  full_time_points <- seq(0, max_follow_up, by = time_interval)
  
  # Generate individual-level parameters
  a0_values <- rnorm(n_individuals, mean = coef_estimates[[1]], sd = coef_se[[1]])
  lambda_i_values <- rnorm(n_individuals, mean = coef_estimates[[2]], sd = coef_se[[2]])
  gamma_values <- runif(n_individuals, min = 0.01, max = 0.05)
  
  # Generate failure times using exponential distribution logic
  fail_times <- log(runif(n_individuals, min = 0, max = 1)) / (-1 * failure_prob)
  failed_flags <- fail_times < max_follow_up
  
  # Preallocate results
  decay_data_list <- vector("list", n_individuals)
  follow_up_years <- rep(NA, n_individuals)
  
  for (i in 1:n_individuals) {
    a0 <- a0_values[i]
    lambda_i <- lambda_i_values[i]
    gamma <- gamma_values[i]
    failure_time <- fail_times[i]
    has_failed <- failed_flags[i]
    
    # Time points depending on failure status
    if (has_failed) {
      censor_time <- failure_time + 3 * time_interval
      times <- full_time_points[full_time_points <= censor_time]
    } else {
      times <- full_time_points
      censor_time <- max_follow_up
    }
    follow_up_years[i] <- max(times)
    
    # Compute true ODn values
    odn_true <- ifelse(
      times < failure_time,
      a0 * exp(-lambda_i * times),
      pmin(a0, a0 * (1 - (1 - exp(-lambda_i * failure_time)) * exp(-gamma * times)))
    )
    
    # Add noise
    noise_sd <- if (noise_model == 1) {
      pmax(sigma_0 + slope_sigma * odn_true, 0.01)
    } else {
      pmax(baseline_noise + fraction * odn_true, 0.01)
    }
    
    noisy_odn <- odn_true
    if (length(times) > 1) {
      noisy_odn[-1] <- pmax(rnorm(length(times) - 1, mean = odn_true[-1], sd = noise_sd[-1]), 0.001)
    }
    
    decay_data_list[[i]] <- data.frame(
      individual = i,
      time = times,
      value = noisy_odn
    )
  }
  
  decay_data <- do.call(rbind, decay_data_list)
  model_parameters <- data.frame(
    id = 1:n_individuals,
    a = a0_values,
    b = lambda_i_values,
    lambda_i = lambda_i_values,
    gamma = gamma_values,
    failure_time = fail_times,
    treatment_failure = failed_flags,
    follow_up_years = follow_up_years
  )
  
  return(list(
    decay_data = decay_data,
    model_parameters = model_parameters
  ))
}


# Call simulate_ODn_3
results <- simulate_ODn_3(
  coef_estimates = coef_estimates,
  coef_se = coef_se,
  n_individuals = 1000,
  sigma_0 = 0.2,
  slope_sigma = 0.1,
  baseline_noise = 0.1,
  fraction = 0.1,
  # failure_prob = 0.2,
  max_follow_up = 10,
  time_interval = 0.5,
  noise_model = 1#,
  # seed = 123
)

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
