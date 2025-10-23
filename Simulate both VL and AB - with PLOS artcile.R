generate_patient_params <- function(
    n = 10000,
    rebound_prop = 0.1,
    min_reb_wk = 26,
    max_reb_wk = 156,
    slope_range_ab_rebound = c(-0.229865, 0.262241),
    slope_range_ab_suppressed = c(-0.0725905, 0.1073314),
    intercept_ab_rebound = c(-10.3950, 5.5760),
    intercept_ab_suppressed = c(-14.5996, 20.4527)
) {
  set.seed(123)
  
  rebound <- rbinom(n, 1, rebound_prop)
  
  params <- data.frame(
    id = 1:n,
    rebound = rebound,
    t_reb = ifelse(rebound == 1,
                   runif(n, min_reb_wk, max_reb_wk),
                   NA),
    # viral slope
    vl_slope = runif(n, 0.05, 0.15),
    # antibody slopes and intercepts by rebound status
    ab_slope = ifelse(rebound == 1,
                      runif(n, slope_range_ab_rebound[1], slope_range_ab_rebound[2]),
                      runif(n, slope_range_ab_suppressed[1], slope_range_ab_suppressed[2])),
    ab_intercept = ifelse(rebound == 1,
                          runif(n, intercept_ab_rebound[1], intercept_ab_rebound[2]),
                          runif(n, intercept_ab_suppressed[1], intercept_ab_suppressed[2])),
    ab_beta = runif(n, 0.01, 0.05),   # antibody rate constant
    ab_lag = runif(n, 4, 12),         # lag in weeks
    ab_baseline = runif(n, 0.5, 1.5)  # baseline antibody level
  )
  
  return(params)
}

simulate_patient_trajectories <- function(params,
                                          times = seq(0, 156, by = 13),  # quarterly (3 years)
                                          detect_threshold = 1000,
                                          max_vl = 1e6,
                                          A_max_fail = 4,
                                          A_max_suppressed = 1.5,
                                          noise_sd_ab = 0.05,
                                          seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  
  sim_data <- vector("list", nrow(params))
  
  for (i in seq_len(nrow(params))) {
    p <- params[i, ]
    
    # --- Viral load trajectory ---
    vl <- sapply(times, function(t) {
      if (p$rebound == 0) {
        runif(1, 1, detect_threshold)  # always suppressed
      } else {
        if (t < p$t_reb) {
          runif(1, 1, detect_threshold / 2)
        } else {
          log10_val <- log10(detect_threshold) + p$vl_slope * (t - p$t_reb)
          10^pmin(log10_val, log10(max_vl))
        }
      }
    })
    
    log_vl <- log10(pmax(vl, 1))
    
    # --- Antibody trajectory (lagged response) ---
    ab <- numeric(length(times))
    ab[1] <- p$ab_baseline
    
    A_max <- ifelse(p$rebound == 1, A_max_fail, A_max_suppressed)
    
    for (t_idx in 2:length(times)) {
      window_start <- max(0, times[t_idx] - p$ab_lag)
      window_idx <- which(times >= window_start & times <= times[t_idx])
      lagged_log_vl <- mean(log_vl[window_idx])
      
      scaled_vl <- pmin(lagged_log_vl / log10(max_vl), 1)
      
      ab[t_idx] <- ab[t_idx - 1] +
        p$ab_beta * (A_max - ab[t_idx - 1]) * scaled_vl +
        rnorm(1, 0, noise_sd_ab)
    }
    
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

# 1. Generate patient parameters
params <- generate_patient_params(n = 1000)

# 2. Simulate viral and antibody trajectories
sim_df <- simulate_patient_trajectories(params)

# 3. Check simulated data
head(sim_df)

# Plot side by side with colors for failing vs suppressed
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

# Antibody plot
p2 <- ggplot(sim_df, aes(x = week, y = antibody, group = id, color = factor(rebound))) +
  geom_line(alpha = 0.4) +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Suppressed", "Failing"),
                     name = "Status") +
  labs(title = "Antibody Trajectories (Patient-Specific Max Levels)",
       x = "Weeks since ART start",
       y = "Antibody level (arbitrary units)") +
  theme_minimal()

# Combine side by side
p1 + p2

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

run_antibody_detection_interval <- function(sim_df, params, 
                                            interval = c("weekly", "biannual", "annual"),
                                            z_threshold = 1.96,
                                            sd_option = "fixed") {
  interval <- match.arg(interval)
  
  # Subset sim_df by chosen interval
  week_seq <- switch(interval,
                     weekly = unique(sim_df$week),
                     biannual = seq(0, max(sim_df$week, na.rm = TRUE), 26),
                     annual = seq(0, max(sim_df$week, na.rm = TRUE), 52))
  sim_sub <- sim_df %>% filter(week %in% week_seq)
  
  # Detect antibody upticks
  antibody_flags <- detect_antibody_upticks(sim_sub,
                                            z_threshold = z_threshold,
                                            sd_option = sd_option)
  
  # Earliest uptick per patient
  ab_summary <- antibody_flags %>%
    filter(uptick_flag == TRUE) %>%
    group_by(id) %>%
    summarise(AB_suspected_failure_time = min(uptick_time, na.rm = TRUE)) %>%
    ungroup()
  
  # Merge with patient parameters
  comparison_df <- params %>%
    left_join(ab_summary, by = "id") %>%
    mutate(
      time_diff = AB_suspected_failure_time - t_reb,
      detected = !is.na(AB_suspected_failure_time),
      true_rebound = !is.na(t_reb),
      correct_detection = case_when(
        detected & true_rebound & AB_suspected_failure_time <= t_reb ~ "FP_early",
        detected & true_rebound & AB_suspected_failure_time > t_reb  ~ "TP",
        detected & !true_rebound                                     ~ "FP",
        !detected & true_rebound                                     ~ "FN",
        TRUE                                                         ~ "TN"
      )
    )
  
  # Summary stats for delays
  summary_stats <- comparison_df %>%
    filter(!is.na(time_diff)) %>%
    summarise(
      mean_delay = mean(time_diff, na.rm = TRUE),
      sd_delay = sd(time_diff, na.rm = TRUE),
      median_delay = median(time_diff, na.rm = TRUE),
      n_events = n()
    )
  
  # 2x2 detection table (using gtsummary)
  detection_table <- comparison_df %>%
    mutate(detected = as.character(detected)) %>%
    select(detected, true_rebound) %>%
    tbl_summary(by = true_rebound)
  
  # Compute diagnostic metrics
  metrics <- comparison_df %>%
    summarise(
      TP = sum(correct_detection == "TP", na.rm = TRUE),
      FP = sum(correct_detection %in% c("FP", "FP_early"), na.rm = TRUE),
      FN = sum(correct_detection == "FN", na.rm = TRUE),
      TN = sum(correct_detection == "TN", na.rm = TRUE)
    ) %>%
    mutate(
      Sensitivity = TP / (TP + FN),
      Specificity = TN / (TN + FP),
      PPV = TP / (TP + FP),
      NPV = TN / (TN + FN),
      Accuracy = (TP + TN) / (TP + TN + FP + FN),
      MeanDelay = summary_stats$mean_delay,
      Scenario = interval,
      LAg_assays_per_year = case_when(
        interval == "weekly"   ~ 52 * nrow(params),
        interval == "biannual" ~ 2  * nrow(params),
        interval == "annual"   ~ 1  * nrow(params)
      ),
      VL_confirmatory_tests = TP + FP
    ) %>%
    select(
      Scenario, LAg_assays_per_year, VL_confirmatory_tests, TP, FP, FN, TN,
      Sensitivity, Specificity, PPV, NPV, Accuracy, MeanDelay
    )
  
  # Return all outputs
  list(
    interval = interval,
    antibody_flags = antibody_flags,
    comparison_df = comparison_df,
    summary_stats = summary_stats,
    detection_table = detection_table,
    metrics = metrics
  )
}

# Example usage:
res_weekly <- run_antibody_detection_interval(sim_df, params, interval = "weekly")
res_biannual <- run_antibody_detection_interval(sim_df, params, interval = "biannual")
res_annual <- run_antibody_detection_interval(sim_df, params, interval = "annual")

# Access outputs:
res_weekly$metrics
res_weekly$detection_table
res_biannual$metrics
res_biannual$detection_table
res_annual$metrics
res_annual$detection_table
