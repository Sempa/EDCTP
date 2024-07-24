# Load necessary libraries
library(dplyr)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# Parameters
n_patients <- 1000               # Number of patients
n_timepoints <- 30               # Number of biannual time points (15 years of data)
baseline_viral_load_mean <- 7    # Mean baseline log viral load
baseline_viral_load_sd <- .01      # Standard deviation of baseline log viral load
decay_rate <- 0.5                # Decay rate for the exponential function
noise_sd <- 0.5                  # Standard deviation of random noise
blip_probability <- 0.2          # Probability of viral rebound between 3.5 to 5 years
blip_magnitude <- 3              # Magnitude of viral blips
decline_probability <- 0.9       # Probability that a patient will decline to <1.3 within 6 months
immediate_decline_probability <- 0.8 # Probability of immediate decline after baseline among those who decline
delayed_decline_probability <- 0.1 # Probability that delayed decliners decline to <1.3 in the next 6 months after one year
occasional_blip_probability <- 0.1  # Probability of occasional viral blips for those who achieve viral load < 1.3
min_viral_load <- 1.3            # Minimum viral load after suppression

# Additional Parameters
mean_age <- 30
sd_age <- 2
female_probability <- 0.6

# Create a data frame to hold patient demographics
patients <- data.frame(
  patient_id = 1:n_patients,
  age = rnorm(n_patients, mean = mean_age, sd = sd_age),
  sex = ifelse(runif(n_patients) < female_probability, "Female", "Male")
)

# Create a data frame to hold the simulated data
sim_data <- data.frame(
  patient_id = rep(1:n_patients, each = n_timepoints),
  timepoint = rep(1:n_timepoints, times = n_patients)
)

# Merge with patient demographics
sim_data <- merge(sim_data, patients, by = "patient_id")

# Simulate viral load data
sim_data <- sim_data %>%
  group_by(patient_id) %>%
  mutate(
    baseline_load = rnorm(1, mean = baseline_viral_load_mean, sd = baseline_viral_load_sd),
    decline_decision = runif(1) < decline_probability,
    immediate_decline = ifelse(decline_decision, runif(1) < immediate_decline_probability, FALSE),
    delay = ifelse(!decline_decision, rbinom(1, 1, delayed_decline_probability), 0),
    exponential_decay = exp(-decay_rate * (timepoint - 1)),
    decay_factor = ifelse(immediate_decline & timepoint == 2, exp(-decay_rate * 1), 
                          ifelse(decline_decision & timepoint <= 2, exp(-decay_rate * (timepoint - 1)), 
                                 ifelse(timepoint > delay * 2 + 2, exp(-decay_rate * (timepoint - (delay * 2 + 2))), 1))),
    noise = rnorm(n(), mean = 0, sd = noise_sd),
    log_viral_load = pmax(baseline_load * exponential_decay * decay_factor + noise, min_viral_load),
    # Gender-specific viral suppression
    log_viral_load = ifelse(sex == "Male" & timepoint <= 4, log_viral_load + 0.5, log_viral_load),
    # Viral rebounds between 3.5 to 5 years
    blip = ifelse(runif(1) < blip_probability & timepoint > 7 & timepoint <= 10, 
                  rnorm(1, mean = blip_magnitude, sd = 0.5), 
                  0),
    log_viral_load = log_viral_load + blip,
    # Suppression after first viral rebound for 15% of those with a rebound
    log_viral_load = ifelse(blip > 0 & runif(1) < 0.15 & timepoint > 10, min_viral_load, log_viral_load),
    # Occasional viral blips
    occasional_blip = ifelse(runif(1) < occasional_blip_probability & log_viral_load < min_viral_load & timepoint > 2, 
                             rnorm(1, mean = blip_magnitude, sd = 0.5), 
                             0),
    log_viral_load = log_viral_load + occasional_blip
  ) %>%
  ungroup()

# Visualize the simulated data
ggplot(sim_data, aes(x = timepoint, y = log_viral_load, group = patient_id, color = sex)) +
  geom_line(alpha = 0.2) +
  labs(title = "Simulated Longitudinal HIV Viral Load Data (Exponential Decay with Blips and Gender Differences)",
       x = "Biannual Time Point",
       y = "Log Viral Load") +
  theme_minimal()
