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
library(patchwork)
library(tidyr)
library(pROC)

model_ab_data <- rbind(readRDS("data/model_data.rds") %>%
                      mutate(phase = "suppressed"),
                    readRDS("data/model_data_from_suppressed.rds")[,c(-6,-7)] %>%
                      mutate(phase = "rebound"))
summary(as.numeric((model_ab_data %>%
           filter(phase == "suppressed"))$intercept_link_log)) ##intercept among suppressed
summary(as.numeric((model_ab_data %>%
           filter(phase == "suppressed"))$slope_link_log) * 7) ##intercept among suppressed
summary(as.numeric((model_ab_data %>%
           filter(phase == "rebound"))$intercept_link_log)) ##intercept among rebound
summary(as.numeric((model_ab_data %>%
           filter(phase == "rebound"))$slope_link_log) * 7) ##intercept among rebound

generate_patient_params <- function(n = 1000,
                                    rebound_prop = 0.10,
                                    min_reb_wk = 0.6,
                                    max_reb_wk = 5,
                                    slope_range_vl = c(0.15, 0.16),
                                    ab_baseline_range = c(0.5, 1.5),
                                    ab_beta_range = c(0.1, 0.5), # AB growth rate
                                    ab_lag_range = c(2, 4)) { # AB growth time lag in comparison to the VL
  ids <- 1:n
  rebound_status <- rbinom(n, 1, rebound_prop)
  
  params <- data.frame(
    id = ids,
    rebound = rebound_status,
    # Viral rebound time (weeks)
    t_reb = ifelse(rebound_status == 1,
                   runif(n, min_reb_wk, max_reb_wk),
                   NA),
    # Viral load slope (log10 scale)
    vl_slope = ifelse(rebound_status == 1,
                      runif(n, slope_range_vl[1], slope_range_vl[2]),
                      runif(n, 0.01, 0.05)),
    # Antibody parameters
    ab_baseline = runif(n, ab_baseline_range[1], ab_baseline_range[2]),
    ab_beta = runif(n, ab_beta_range[1], ab_beta_range[2]),
    ab_lag = runif(n, ab_lag_range[1], ab_lag_range[2])
  )
  
  return(params)
}

### 2. Simulation function with separate max antibody levels
simulate_patient_trajectories <- function(params,
                                          times = seq(0, 156, by = 13),  # quarterly for 3 years
                                          detect_threshold = 1e3,
                                          max_vl = 1e6,
                                          A_max_fail = 4,          # max for failing patients
                                          A_max_suppressed = 1.5,  # max for suppressed patients
                                          noise_sd_ab = 0.05,
                                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sim_data <- vector("list", nrow(params))
  
  for (i in seq_len(nrow(params))) {
    p <- params[i, ]
    
    # --- Viral load trajectory ---
    vl <- sapply(times, function(t) {
      if (p$rebound == 0) {
        runif(1, 1, detect_threshold)  # always suppressed
      } else {
        if (t < p$t_reb * 7) {
          runif(1, 1, detect_threshold / 2)
        } else {
          log10_val <- log10(detect_threshold) +
            p$vl_slope * (t - p$t_reb * 7)
          10^pmin(log10_val, log10(max_vl))
        }
      }
    })
    
    log_vl <- log10(pmax(vl, 1))
    
    # --- Antibody trajectory: saturating response ---
    ab <- numeric(length(times))
    ab[1] <- p$ab_baseline  # starting antibody value
    
    # choose max antibody per patient
    A_max <- ifelse(p$rebound == 1, A_max_fail, A_max_suppressed)
    
    for (t_idx in 2:length(times)) {
      # lagged log VL
      window_start <- max(0, times[t_idx] - p$ab_lag)
      window_idx <- which(times >= window_start & times <= times[t_idx])
      lagged_log_vl <- mean(log_vl[window_idx])
      
      # scale lagged VL to 0-1 for gradual growth
      scaled_vl <- pmin(lagged_log_vl / log10(max_vl), 1)
      
      # logistic/saturating growth toward patient-specific A_max
      ab[t_idx] <- ab[t_idx - 1] + p$ab_beta * (A_max - ab[t_idx - 1]) * scaled_vl +
        rnorm(1, 0, noise_sd_ab)
    }
    
    # constrain antibody within 0 and patient-specific max
    ab <- pmin(pmax(ab, 0.001), A_max)
    
    sim_data[[i]] <- data.frame(
      id = p$id,
      week = times,
      vl = vl,
      log_vl = log_vl,
      antibody = ab,
      rebound = p$rebound
    )
  }
  
  sim_df <- do.call(rbind, sim_data)
  return(sim_df)
}

### 3. Run simulation
params <- generate_patient_params(n = 1000, rebound_prop = 0.15)
sim_df <- simulate_patient_trajectories(params, seed = 42)

### 4. Visualization
# Viral Load
p1 <- ggplot(sim_df, aes(x = week, y = vl, group = id, color = factor(rebound))) +
  geom_line(alpha = 0.3) +
  geom_hline(yintercept = 1000, linetype = "dashed", color = "black") +
  scale_y_log10() +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Suppressed", "Failing"),
                     name = "Status") +
  labs(title = "HIV Viral Load Trajectories",
       x = "Weeks since ART start",
       y = "Viral load (copies/mL, log10 scale)") +
  theme_minimal()

# Antibody
p2 <- ggplot(sim_df, aes(x = week, y = antibody, group = id, color = factor(rebound))) +
  geom_line(alpha = 0.4) +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Suppressed", "Failing"),
                     name = "Status") +
  labs(title = "Antibody Trajectories (Patient-Specific Max Levels)",
       x = "Weeks since ART start",
       y = "Antibody level (arbitrary units)") +
  theme_minimal()

p1 + p2

# # Detect antibody upticks
# detect_antibody_upticks <- function(
#     antibody_data,
#     z_threshold = 1.96,
#     min_history = 2,          # applies only if sd_option = "expanding"
#     sd_option = c("fixed", "expanding"), 
#     fixed_sd = 0.3141,        # assay variability
#     seed = 123
# ) {
#   set.seed(seed)
#   sd_option <- match.arg(sd_option)
#   
#   # adjust min_history dynamically
#   min_history_used <- if (sd_option == "fixed") 1 else min_history
#   
#   antibody_data <- antibody_data %>%
#     arrange(id, week) %>%
#     group_by(id) %>%
#     mutate(
#       z_score = sapply(seq_along(antibody), function(i) {
#         if (i <= min_history_used) return(NA_real_)
#         
#         # use previous values to compute the historical mean
#         past_vals <- antibody[1:(i - 1)]
#         mean_past <- mean(past_vals, na.rm = TRUE)
#         
#         # define SD according to chosen method
#         sd_val <- switch(sd_option,
#                          fixed = fixed_sd,
#                          expanding = sd(past_vals, na.rm = TRUE))
#         
#         if (is.na(sd_val) || sd_val == 0) return(NA_real_)
#         
#         # compute z-score
#         (antibody[i] - mean_past) / sd_val
#       }),
#       
#       # flag significant upticks
#       uptick_flag = z_score > z_threshold
#     ) %>%
#     mutate(
#       # identify the first uptick per patient
#       uptick_time = ifelse(row_number() == min(which(uptick_flag), na.rm = TRUE), week, NA_real_)
#     ) %>%
#     fill(uptick_time, .direction = "down") %>%
#     ungroup()
#   
#   return(antibody_data)
# }

detect_antibody_upticks <- function(
    antibody_data,
    z_threshold = 1.96,
    min_history = 2,                   # used only for expanding SD
    sd_option = c("fixed", "expanding", "heteroskedastic"),
    fixed_sd = 0.3141,                 # for fixed SD
    seed = 123
) {
  set.seed(seed)
  sd_option <- match.arg(sd_option)
  
  # Adjust minimum history requirement based on SD method
  min_history_used <- if (sd_option %in% c("fixed", "heteroskedastic")) 1 else min_history
  
  antibody_data <- antibody_data %>%
    arrange(id, week) %>%
    group_by(id) %>%
    mutate(
      z_score = sapply(seq_along(antibody), function(i) {
        if (i <= min_history_used) return(NA_real_)
        
        # Past antibody values up to the current timepoint
        past_vals <- antibody[1:(i - 1)]
        mean_past <- mean(past_vals, na.rm = TRUE)
        
        # Choose SD calculation method
        sd_val <- switch(
          sd_option,
          fixed = fixed_sd,
          expanding = sd(past_vals, na.rm = TRUE),
          heteroskedastic = max(0.001, -0.01469 + 0.14513 * mean_past)  # Linear SD model
        )
        
        if (is.na(sd_val) || sd_val <= 0) return(NA_real_)
        
        # Compute Z-score: how many SDs above the mean is the current value?
        (antibody[i] - mean_past) / sd_val
      }),
      
      # Flag upticks that exceed the threshold
      uptick_flag = z_score > z_threshold
    ) %>%
    mutate(
      # Mark the first time a significant uptick occurs per patient
      uptick_time = ifelse(
        row_number() == min(which(uptick_flag), na.rm = TRUE),
        week,
        NA_real_
      )
    ) %>%
    fill(uptick_time, .direction = "down") %>%
    ungroup()
  
  return(antibody_data)
}

# Apply to simulated data
antibody_flags <- detect_antibody_upticks(sim_df, 
                                          z_threshold = 1.96, 
                                          sd_option = "fixed") #heteroskedastic,expanding


six_monthly_data <- detect_antibody_upticks(sim_df %>%
                                              filter(week %in% seq(0,156, 26)), 
                                            z_threshold = 1.96, 
                                            sd_option = "fixed")
# Inspect flagged individuals
head(antibody_flags %>% filter(uptick_flag == TRUE))

# Summarize: earliest uptick per patient
ab_summary <- antibody_flags %>%
  filter(uptick_flag == TRUE) %>%
  group_by(id) %>%
  summarise(AB_suspected_failure_time = min(uptick_time, na.rm = TRUE)) %>%
  ungroup()

# Merge with patient parameters (t_reb = viral rebound time)
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

# Quick summary of ordering
table(comparison_df$event_order, useNA = "ifany")

# Summary stats for time lag (patients with both events)
comparison_df %>%
  filter(!is.na(time_diff)) %>%
  summarise(
    mean_diff = mean(time_diff, na.rm = TRUE),
    sd_diff = sd(time_diff, na.rm = TRUE),
    median_diff = median(time_diff, na.rm = TRUE),
    n = n()
  )

# Visualization: antibody vs viral rebound timing
ggplot(comparison_df, aes(x = t_reb, y = AB_suspected_failure_time, color = event_order)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x = "Time of Viral Rebound (t_reb)",
    y = "Time of Suspected Antibody Uptick",
    title = "Comparison of Viral Rebound and Antibody Uptick Timing",
    subtitle = "Dashed line = same timing; Points above line = antibody uptick occurs later",
    color = "Event Order"
  ) +
  theme_minimal(base_size = 13)

# Derive detection classification
comparison_df <- comparison_df %>%
  mutate(
    detected = !is.na(AB_suspected_failure_time),
    true_rebound = !is.na(t_reb),
    correct_detection = case_when(
      detected & true_rebound & AB_suspected_failure_time <= t_reb ~ "TP",
      detected & true_rebound & AB_suspected_failure_time > t_reb ~ "FP_late",
      detected & !true_rebound ~ "FP",
      !detected & true_rebound ~ "FN",
      TRUE ~ "TN"
    )
  )

# Tabulate results
table_results <- table(comparison_df$correct_detection)

# Summarise performance metrics
performance <- tibble(
  TP = sum(comparison_df$correct_detection == "TP"),
  FP = sum(comparison_df$correct_detection %in% c("FP", "FP_late")),
  FN = sum(comparison_df$correct_detection == "FN"),
  TN = sum(comparison_df$correct_detection == "TN")
) %>%
  mutate(
    Sensitivity = TP / (TP + FN),
    Specificity = TN / (TN + FP),
    PPV = TP / (TP + FP),
    NPV = TN / (TN + FN),
    Accuracy = (TP + TN) / (TP + FP + FN + TN)
  )

performance

ab_summary2 <- antibody_flags %>%
  group_by(id) %>%
  summarise(max_antibody = max(antibody, na.rm = TRUE)) %>%
  ungroup()

roc_data <- ab_summary2 %>%
  left_join(params %>% select(id, t_reb, rebound), by = "id") %>%
  mutate(
    rebound_flag = rebound == 1,   # TRUE if rebound
    valid = !is.na(max_antibody) & !is.na(rebound_flag)
  ) %>%
  filter(valid)

roc_obj <- roc(
  response = roc_data$rebound_flag,
  predictor = roc_data$max_antibody,
  direction = "<"  # lower antibody → less likely rebound
)

plot(roc_obj, col = "blue", main = "ROC: Antibody vs Viral Rebound")
auc(roc_obj)

x <- comparison_df %>%
  mutate(detected = as.character(detected)) %>%
  dplyr::select(detected, true_rebound) %>%
  tbl_summary(by = true_rebound)
x
y <- comparison_df %>%
  filter(!is.na(t_reb) & !is.na(AB_suspected_failure_time))
mean(y$AB_suspected_failure_time)
sd(y$AB_suspected_failure_time)

run_antibody_detection <- function(sim_df, params,
                                   interval = c("weekly", "biannual", "annual"),
                                   z_threshold = 1.96,
                                   sd_option = "fixed") {
  library(dplyr)
  library(gtsummary)
  
  # Choose interval
  interval <- match.arg(interval)
  week_seq <- switch(interval,
                     weekly = unique(sim_df$week),
                     biannual = seq(0, max(sim_df$week, na.rm = TRUE), 26),
                     annual = seq(0, max(sim_df$week, na.rm = TRUE), 52))
  
  # Subset based on chosen interval
  sim_sub <- sim_df %>%
    filter(week %in% week_seq)
  
  # Detect upticks
  antibody_flags <- detect_antibody_upticks(sim_sub,
                                            z_threshold = z_threshold,
                                            sd_option = sd_option)
  
  # Summarize earliest uptick per patient
  ab_summary <- antibody_flags %>%
    filter(uptick_flag == TRUE) %>%
    group_by(id) %>%
    summarise(AB_suspected_failure_time = min(uptick_time, na.rm = TRUE)) %>%
    ungroup()
  
  # Merge with patient parameters (t_reb = viral rebound time)
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
  
  # Summary stats
  summary_stats <- comparison_df %>%
    filter(!is.na(time_diff)) %>%
    summarise(
      mean_diff = mean(time_diff, na.rm = TRUE),
      sd_diff = sd(time_diff, na.rm = TRUE),
      median_diff = median(time_diff, na.rm = TRUE),
      n = n()
    )
  
  # Derive detection classification
  comparison_df <- comparison_df %>%
    mutate(
      detected = !is.na(AB_suspected_failure_time),
      true_rebound = !is.na(t_reb),
      correct_detection = case_when(
        detected & true_rebound & AB_suspected_failure_time <= t_reb ~ "TP",
        detected & true_rebound & AB_suspected_failure_time > t_reb ~ "FP_late",
        detected & !true_rebound ~ "FP",
        !detected & true_rebound ~ "FN",
        TRUE ~ "TN"
      )
    )
  
  # Tabulate results
  table_results <- table(comparison_df$correct_detection)
  
  # Create summary table
  gtsummary_tbl <- comparison_df %>%
    mutate(detected = as.character(detected)) %>%
    select(detected, true_rebound) %>%
    tbl_summary(by = true_rebound)
  
  # Return all results
  list(
    antibody_flags = antibody_flags,
    comparison_df = comparison_df,
    summary_stats = summary_stats,
    table_results = table_results,
    gtsummary = gtsummary_tbl
  )
}

# Weekly analysis
res_weekly <- run_antibody_detection(sim_df, params, interval = "weekly")

# Biannual (6-month) analysis
res_biannual <- run_antibody_detection(sim_df, params, interval = "biannual")

# Annual analysis
res_annual <- run_antibody_detection(sim_df, params, interval = "annual")

# Access results
## weekly
res_weekly$summary_stats
res_weekly$table_results
res_weekly$gtsummary
## Biannual
res_biannual$summary_stats
res_biannual$table_results
res_biannual$gtsummary
## Annual
res_annual$summary_stats
res_annual$table_results
res_annual$gtsummary