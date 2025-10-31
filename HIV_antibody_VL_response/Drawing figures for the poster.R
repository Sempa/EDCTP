library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# --- Parameter grid ---
spec_values <- c(0.90, 0.95, 0.97)
AB_rebound_delay_seq <- seq(0, 0.4, by = 0.01)
annual_rebound_rate_seq <- seq(0.04, 0.5, by = 0.01)

param_grid <- expand.grid(
  AB_Specificity = spec_values,
  AB_rebound_delay = AB_rebound_delay_seq,
  annual_rebound_rate = annual_rebound_rate_seq
)

population <- 1000

# --- Function that computes metrics for one parameter set ---
compute_metrics <- function(AB_Specificity, AB_rebound_delay, annual_rebound_rate, population = 1000) {
  monitoring_scenario <- c('annual VL', 'Annual AB', 'biannual AB')
  N <- rep(population, length(monitoring_scenario))
  prim_VL_tests <- c(population, 0, 0)
  testing_interval <- c(1, 1, 0.5)
  AB_testing <- c(0, rep(population, 2))
  
  number_of_rebounds_PerAnnum <- rep(population * annual_rebound_rate, 3)
  number_of_rebounds_Per_test_round <- c(
    NA,
    population * annual_rebound_rate * (testing_interval[2] + AB_rebound_delay),
    (testing_interval[3] * annual_rebound_rate * population) + (population * annual_rebound_rate * AB_rebound_delay)
  )
  
  true_positives <- c(
    population * annual_rebound_rate,
    number_of_rebounds_PerAnnum[2] * testing_interval[2],
    number_of_rebounds_PerAnnum[3] * testing_interval[3]
  )
  
  false_positives <- c(
    0,
    (N[2] - number_of_rebounds_PerAnnum[2]) * (1 - AB_Specificity),
    (N[3] - number_of_rebounds_PerAnnum[3]) * (1 - AB_Specificity)
  )
  
  vl_confirmation <- c(
    0,
    true_positives[2] + false_positives[2],
    true_positives[3] + false_positives[3]
  )
  
  total_VL <- c(
    vl_confirmation[1] + prim_VL_tests[1],
    vl_confirmation[2] + prim_VL_tests[2],
    vl_confirmation[3] + prim_VL_tests[3]
  )
  
  mean_delay_yrs <- c(
    0.5,
    0.5 * (AB_rebound_delay + testing_interval[2] + AB_rebound_delay),
    0.5 * (AB_rebound_delay + testing_interval[3] + AB_rebound_delay)
  )
  mean_delay_months <- mean_delay_yrs * 12
  
  total_VL_per_yr <- c(
    N[1],
    total_VL[2] / testing_interval[2],
    total_VL[3] / testing_interval[3]
  )
  
  tibble(
    Scenario = monitoring_scenario,
    AB_Specificity = AB_Specificity,
    AB_rebound_delay = AB_rebound_delay,
    annual_rebound_rate = annual_rebound_rate,
    Mean_delay_months = mean_delay_months,
    Total_VL_per_year = total_VL_per_yr
  )
}

# --- Run across grid ---
results <- param_grid %>%
  mutate(data = pmap(
    list(AB_Specificity, AB_rebound_delay, annual_rebound_rate),
    ~ compute_metrics(..1, ..2, ..3, population)
  )) %>%
  unnest(data, names_sep = "_")

# --- Plot distributions ---

## 1. Distribution of Mean Delay (months)
ggplot(results, aes(x = data_Mean_delay_months, fill = Scenario)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ data_AB_Specificity, labeller = label_both) +
  labs(
    title = "Distribution of Mean Delay (months)",
    subtitle = "Varying AB rebound delay (0–0.4) and annual rebound rate (0.04–0.5)",
    x = "Mean Delay (months)",
    y = "Density"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

## 2. Distribution of Total VL (per year)
ggplot(results, aes(x = data_Total_VL_per_year, fill = data_Scenario)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ AB_Specificity, labeller = label_both, scales = "free_x") +
  labs(
    title = "Distribution of Total VL (per year)",
    subtitle = "Varying AB rebound delay (0–0.4) and annual rebound rate (0.04–0.5)",
    x = "Total VL (per year)",
    y = "Density"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")
