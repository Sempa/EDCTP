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

model_ab_data <- rbind(readRDS("data/model_data.rds") %>%
                      mutate(phase = "suppressed"),
                    readRDS("data/model_data_from_suppressed.rds")[,c(-6,-7)] %>%
                      mutate(phase = "rebound"))

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

x <- filtered_subjects %>%
  group_by(subject_label_blinded) %>%
  # keep only visits at or before 1 year
  filter(years_since_tx_start <= 2) %>%
  # pick the last visit for each subject
  slice_max(order_by = years_since_tx_start, n = 1, with_ties = FALSE) %>%
  ungroup()
summary(x$sedia_ODn)

simulate_vl_simple <- function(n = 84,
                               times = seq(0, 156, by = 0.1), # weeks (0 to 10 weeks)
                               min_reb_wk = 4/7,  # ≈0.6 weeks
                               max_reb_wk = 35/7, # ≈5 weeks
                               detect_threshold = 200,
                               max_vl = 1e6,
                               slope_range = c(0.2, 0.6), # slope per week on log10 scale
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # random rebound times (in weeks)
  t_reb <- runif(n, min_reb_wk, max_reb_wk)
  # random slopes for linear rise on log10 scale
  slopes <- runif(n, slope_range[1], slope_range[2])
  
  out <- lapply(1:n, function(i) {
    tr <- t_reb[i]
    sl <- slopes[i]
    
    vl <- sapply(times, function(t) {
      if (t < tr) {
        runif(1, 1, detect_threshold/2) # below detection
      } else {
        # log10-linear increase starting from threshold at rebound
        log10_val <- log10(detect_threshold) + sl * (t - tr)
        10^pmin(log10_val, log10(max_vl)) # cap at max_vl
      }
    })
    
    data.frame(id = i, time_wk = times, vl = vl)
  })
  
  do.call(rbind, out)
}

# Example usage
sim <- simulate_vl_simple(seed = 123)
library(ggplot2)
ggplot2::ggplot(data = sim, aes(x = time_wk, y = vl, group = id)) +
  geom_line(alpha = 0.5) +
  # theme_minimal() +
  labs(x = "Weeks since ART stop", y = "Viral load (copies/mL, log10 scale)") +
  theme_minimal()

library(ggplot2)
library(dplyr)

### 1️⃣ Function to generate patient-level parameters ----
generate_patient_params <- function(n = 1000,
                                    rebound_prop = 0.10,
                                    min_reb_wk = 0.6,
                                    max_reb_wk = 5,
                                    slope_range_vl = c(0.2, 0.6),  # VL slope (log10 scale)
                                    slope_range_ab = c(-4.691e-03, 5.352e-03), # Antibody slope (log scale)
                                    intercept_ab = c(-1.48500, 0.79657)) {     # Antibody baseline
  
  ids <- 1:n
  rebound_status <- rbinom(n, 1, rebound_prop)
  
  params <- data.frame(
    id = ids,
    rebound = rebound_status,
    t_reb = ifelse(rebound_status == 1,
                   runif(n, min_reb_wk, max_reb_wk),
                   NA),
    vl_slope = ifelse(rebound_status == 1,
                      runif(n, slope_range_vl[1], slope_range_vl[2]),
                      runif(n, 0.01, 0.05)),  # flatter if suppressed
    ab_slope = ifelse(rebound_status == 1,
                      runif(n, slope_range_ab[1], slope_range_ab[2]),
                      runif(n, -0.05, 0.02)), # decline or flat if suppressed
    ab_intercept = runif(n, intercept_ab[1], intercept_ab[2])
  )
  
  return(params)
}

### 2️⃣ Function to simulate trajectories ----
simulate_patient_trajectories <- function(params,
                                          times = seq(0, 156, by = 13),  # quarterly over 3 years (weeks)
                                          detect_threshold = 200,
                                          max_vl = 1e6,
                                          ab_max = 4,
                                          noise_sd_ab = 0.05,
                                          seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  sim_data <- list()
  
  for (i in 1:nrow(params)) {
    p <- params[i, ]
    
    vl <- sapply(times, function(t) {
      if (p$rebound == 0) {
        runif(1, 1, detect_threshold)  # always suppressed
      } else {
        if (t < p$t_reb * 7) {         # convert rebound time to days
          runif(1, 1, detect_threshold / 2)
        } else {
          log10_val <- log10(detect_threshold) + p$vl_slope * ((t - p$t_reb * 7) / 7)
          10^pmin(log10_val, log10(max_vl))
        }
      }
    })
    
    # antibody trajectories (exponential/log-link)
    log_ab <- p$ab_intercept + p$ab_slope * (times / 13) + rnorm(length(times), 0, noise_sd_ab)
    ab <- exp(log_ab)
    ab <- pmin(ab, ab_max)
    
    sim_data[[i]] <- data.frame(
      id = p$id,
      week = times,
      vl = vl,
      antibody = ab,
      rebound = p$rebound
    )
  }
  
  sim_df <- do.call(rbind, sim_data)
  return(sim_df)
}

### 3️⃣ Run simulation ----
params <- generate_patient_params(n = 1000, rebound_prop = 0.1)
sim_df <- simulate_patient_trajectories(params, seed = 42)

### 4️⃣ Visualization ----

# Viral Load Plot
p1 <- ggplot(sim_df, aes(x = week, y = vl, group = id, color = factor(rebound))) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = 1000, linetype = "dashed", color = "black") +
  scale_y_log10() +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Suppressed", "Rebound"),
                     name = "Status") +
  labs(title = "HIV Viral Load Trajectories",
       x = "Weeks since ART interruption",
       y = "Viral load (copies/mL, log10 scale)") +
  theme_minimal()

# Antibody Plot
p2 <- ggplot(sim_df, aes(x = week, y = antibody, group = id, color = factor(rebound))) +
  geom_line(alpha = 0.3) +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Suppressed", "Rebound"),
                     name = "Status") +
  labs(title = "HIV Antibody Trajectories",
       x = "Weeks since ART interruption",
       y = "Antibody level (arbitrary units)") +
  theme_minimal()

# Display side-by-side
library(patchwork)
p1 + p2

detect_antibody_upticks <- function(
    antibody_data,
    z_threshold = 1.96,
    min_history = 2,        # minimum number of past values required
    sd_option = c("fixed", "expanding"), 
    fixed_sd = 0.15,        # optional fixed SD, can tune based on assay variability
    seed = 123
) {
  set.seed(seed)
  library(dplyr)
  library(tidyr)
  
  sd_option <- match.arg(sd_option)
  
  antibody_data <- antibody_data %>%
    arrange(id, week) %>%
    group_by(id) %>%
    mutate(
      z_score = sapply(seq_along(antibody), function(i) {
        # not enough history yet
        if (i <= min_history) return(NA_real_)
        
        # historical antibody values up to (i-1)
        past_vals <- antibody[1:(i - 1)]
        mean_past <- mean(past_vals, na.rm = TRUE)
        
        # define SD based on chosen option
        sd_val <- switch(sd_option,
                         fixed = fixed_sd,
                         expanding = sd(past_vals, na.rm = TRUE))
        
        if (is.na(sd_val) || sd_val == 0) return(NA_real_)
        
        # compute z-score: how many SDs above the past mean is current value?
        (antibody[i] - mean_past) / sd_val
      }),
      
      uptick_flag = z_score > z_threshold
    ) %>%
    mutate(
      # identify first uptick time
      uptick_time = ifelse(row_number() == which(uptick_flag)[1], week, NA_real_)
    ) %>%
    fill(uptick_time, .direction = "down") %>%
    mutate(
      uptick_time = ifelse(row_number() == 1, NA_real_, uptick_time)
    ) %>%
    ungroup()
  
  return(antibody_data)
}

# Example simulated dataset (from previous antibody simulation)
set.seed(1)
sim_example <- simulate_antibody_simple(n = 10, rebound_data = rebound, phase = "rebound")

# Detect antibody upticks
antibody_flags <- detect_antibody_upticks(sim_df, 
                                          z_threshold = 1.96, 
                                          sd_option = "expanding")

# Inspect flagged individuals
head(antibody_flags %>% filter(uptick_flag == TRUE))

# Summarize antibody flags: get the first (earliest) uptick_time per patient
ab_summary <- antibody_flags %>%
  filter(uptick_flag == TRUE) %>%
  group_by(id) %>%
  summarise(AB_suspected_failure_time = min(uptick_time, na.rm = TRUE)) %>%
  ungroup()

# Merge with params (which has t_reb)
comparison_df <- params %>%
  left_join(ab_summary, by = "id") %>%
  mutate(
    time_diff = AB_suspected_failure_time - t_reb,
    event_order = case_when(
      is.na(AB_suspected_failure_time) ~ "No antibody failure detected",
      AB_suspected_failure_time > t_reb ~ "Rebound before antibody uptick",
      AB_suspected_failure_time < t_reb ~ "Antibody uptick before rebound",
      AB_suspected_failure_time == t_reb ~ "Simultaneous"
    )
  )

# Quick summary of the ordering of events
table(comparison_df$event_order, useNA = "ifany")

# Summary stats for time lag (only those with both events)
comparison_df %>%
  filter(!is.na(time_diff)) %>%
  summarise(
    mean_diff = mean(time_diff, na.rm = TRUE),
    sd_diff = sd(time_diff, na.rm = TRUE),
    median_diff = median(time_diff, na.rm = TRUE),
    n = n()
  )

# Visualization
ggplot(comparison_df, aes(x = t_reb, y = AB_suspected_failure_time, color = event_order)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x = "Time of Viral Rebound (t_reb)",
    y = "Time of Suspected Antibody Failure (AB_suspected_failure_time)",
    title = "Comparison of Viral Rebound and Antibody Uptick Timing",
    subtitle = "Dashed line = same timing; Points above line = antibody failure occurs later",
    color = "Event Order"
  ) +
  theme_minimal(base_size = 13)
