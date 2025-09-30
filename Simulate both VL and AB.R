x <- filtered_subjects %>%
  group_by(subject_label_blinded) %>%
  # keep only visits at or before 1 year
  filter(years_since_tx_start <= 2) %>%
  # pick the last visit for each subject
  slice_max(order_by = years_since_tx_start, n = 1, with_ties = FALSE) %>%
  ungroup()
summary(x$sedia_ODn)

simulate_vl_simple <- function(n = 84,
                               times = seq(0, 10, by = 0.1), # weeks (0 to 10 weeks)
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
