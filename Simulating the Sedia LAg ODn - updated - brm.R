# Load necessary libraries
library(brms)
library(survival)
library(ggplot2)
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
         Sex = BiologicalSex,
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

# Fit the linear mixed effects model
fit <- brm(
  sedia_ODn ~ poly(years_since_tx_start, 2) + (1 + years_since_tx_start | subject_label_blinded),
  data = full_dataset  %>%
    filter(Group == 'early suppressions'),
  family = gaussian(),
  prior = prior(normal(0, 1), class = "b"),
  chains = 4,
  cores = 4,
  warmup = 1000,
  iter = 2000
)

# Extract the fixed effects coefficients
fixed_effects <- fixef(fit)

# Extract the individual-specific random effects
random_effects <- ranef(fit)

# Extract the standard deviation of the random effects
sd_random_effects1 <- summary(fit)$random$subject_label_blinded[2,1]
sd_random_effects2 <- summary(fit)$random$subject_label_blinded[3,1]
# Define the Weibull decay distribution parameters
# Simulate Weibull decay distribution data using the model parameters
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

# Define a two-degree polynomial to generate decay parameters
generate_decay_rate <- function(t) {
  # browser()
  # Two-degree polynomial: a*t^2 + b*t + c
  a <- rnorm(1, mean = fixed_effects[3,1], 
             sd = summary(fit)$random$subject_label_blinded[3,1])#1.818730
  b <- rnorm(1, mean = fixed_effects[2,1], 
             sd = summary(fit)$random$subject_label_blinded[2,1])#-8.97018
  c <- rnorm(1, mean = fixed_effects[1,1], 
             sd = summary(fit)$random$subject_label_blinded[1,1])#3.254208
  decay_rate <- a * t^2 + b * t + c
  return(c(pmax(decay_rate, 0.05), a, b, c))  # Ensure decay rates are non-negative
}

# Generate decay rates for each individual
decay_rates <- sapply(1:n_individuals, function(i) {
  generate_decay_rate(time_points)  # Each individual has a unique decay rate
})

# Define the exponential decay function
exp_decay <- function(t, baseline, decay_rate) {
  pmax(baseline * exp(-decay_rate * t), .001)  # Ensure positive ODn values
}

# Generate decay data for each individual
decay_data <- data.frame(time = rep(time_points, n_individuals),
                         individual = rep(1:n_individuals, each = length(time_points)),
                         value = unlist(lapply(1:n_individuals, function(i) {
                           exp_decay(time_points, baselines[i], decay_rates[1:length(n_individuals), i])
                         })))
model_parameters <- bind_cols(id = 1:n_individuals,
                              baseline = as.data.frame(baselines), 
                              remove_rownames(data.frame(t(data.frame(decay_rates[c(22,23,24),])))) %>%
                                dplyr::select(a = X1, b = X2, c = X3))
# model_parameters <- readRDS('model_parameters_100K.rds')
# decay_data <- readRDS('decay_data_100K.rds')
# decay_rates <- readRDS('decay_rates_100K.rds')
# Plot the decay curves
ggplot(decay_data, aes(x = time, y = value, group = individual)) +
  geom_line(alpha = 0.5) +
  labs(title = "Exponential Decay Curves for Individuals",
       x = "Time since ART start (years)",
       y = "ODn Value") +
  theme_minimal()
