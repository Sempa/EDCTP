# ce_ab_triage.R
# Cost-effectiveness comparing Ab-triage -> PCR vs PCR-only
# Requirements: tidyverse, data.table, ggplot2, cowplot, reshape2
# Install missing packages if needed:
# install.packages(c("tidyverse","data.table","ggplot2","cowplot","reshape2","scales","triangle"))

library(tidyverse)
library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
# triangle package optional for triangular draws
# install.packages("triangle"); library(triangle)

set.seed(12345)

# ----------------------------
# 1) Model inputs (base case)
# ----------------------------
params <- list(
  # Population & epidemic
  N = 10000,                       # cohort size (for scaling)
  annual_non_suppression = 0.10,   # baseline annual proportion with VL non-suppression (viremic) - user varies this
  
  # Test characteristics (base case)
  ab_sensitivity = 0.6,            # sensitivity of Ab test to detect recent VL rebound (fraction)
  ab_specificity = 0.95,           # specificity
  pcr_sensitivity = 1.0,           # assumed
  pcr_specificity = 1.0,           # assumed
  
  # Costs (USD)
  cost_ab = 2.0,                   # cost per Ab test
  cost_pcr = 40.0,                 # cost per PCR test
  cost_clinic_visit = 5.0,         # incremental cost if visit needed for test/triage (per test round)
  
  # Effectiveness: benefit of earlier detection
  # For a person whose VL rebound is detected earlier via Ab triage, we model an expected gain (years of VL suppression regained earlier)
  # This is a shorthand: you should replace with a dynamic model if you have longitudinal trajectories.
  mean_effect_years_gain_on_detection = 0.25,  # e.g. 3 months earlier suppression on average per detected rebound
  # stratify resistance vs behavioral non-adherence:
  prop_reistance = 0.4,            # proportion of viremic who have resistance (vs behaviour)
  p_resupp_after_resistance_detect = 0.2,   # probability of re-suppression if resistance detected (may require regimen change)
  p_resupp_after_behav_detect = 0.6,       # probability of re-suppression if behavioural non-adherence detected
  # multiplier for earlier detection efficacy:
  efficacy_multiplier_resistance = 0.6,    # detection less effective for resistance (lower chance to re-suppress)
  efficacy_multiplier_behav = 1.0          # detection more effective for behaviour
)

# Example WTP thresholds (user originally requested $200, $500, $800)
WTPs <- c(200, 500, 800)

# ----------------------------
# 2) Helper functions
# ----------------------------

# Expected annual costs and effectiveness for a single person under given monitoring strategy
# frequency_yrs: Ab testing interval in years (e.g., 1 = annual, 0.5 = biannual, 0.0833 ~ monthly)
# strategy: "PCR_only" or "AB_triage"
calc_per_person <- function(frequency_yrs = 1,
                            params_list = params,
                            annual_non_suppression = params_list$annual_non_suppression,
                            ab_sens = params_list$ab_sensitivity,
                            ab_spec = params_list$ab_specificity,
                            cost_ab = params_list$cost_ab,
                            cost_pcr = params_list$cost_pcr,
                            cost_visit = params_list$cost_clinic_visit,
                            mean_gain = params_list$mean_effect_years_gain_on_detection,
                            prop_res = params_list$prop_reistance,
                            p_resup_res = params_list$p_resupp_after_resistance_detect,
                            p_resup_beh = params_list$p_resupp_after_behav_detect,
                            eff_mult_res = params_list$efficacy_multiplier_resistance,
                            eff_mult_beh = params_list$efficacy_multiplier_behav
) {
  
  # expected number of Ab tests per year
  ab_tests_per_year <- ifelse(is.infinite(frequency_yrs) || frequency_yrs <= 0, 0, 1 / frequency_yrs)
  
  # PCR-only: assume 1 PCR per person per year (or 0? choose policy)
  # We'll assume PCR-only means routine PCR once per year for everyone (you can change).
  pcr_tests_per_year_pcronly <- 1   # changeable if needed
  
  # Costs under PCR-only per person-year
  cost_pcronly <- pcr_tests_per_year_pcronly * (cost_pcr + cost_visit)
  
  # Effectiveness: PCR-only assumed to detect all rebounds at that routine PCR visit; average delay to detection = 0.5 year if annual
  # For simplicity, assume average delay under PCR-only = testing_interval/2 = 0.5 year
  # But many rebounds may happen between visits; for comparability we compute expected "delay" to detection:
  delay_pcronly_years <- 0.5 * (1 / pcr_tests_per_year_pcronly)
  # expected gain in suppression-year if detected earlier compared to no detection:
  # we treat the gain as mean_gain * probability of resuppression on detection (weighted)
  p_resupp_weighted <- prop_res * p_resup_res + (1 - prop_res) * p_resup_beh
  eff_pcronly <- annual_non_suppression * p_resupp_weighted * mean_gain  # years-of-suppression gained per person-year
  
  # AB_triage:
  # If Ab test positive -> confirm with PCR (assume PCR applied only to positives)
  # Expected numbers:
  ab_pos_if_viremic <- ab_sens       # sensitivity among viremic
  ab_pos_if_not <- 1 - ab_spec       # false positive rate among non-viremic
  
  # expected PCR confirmations per person-year under triage:
  # P(ab positive) = P(viremic)*ab_sens + (1-P(viremic))*(1-ab_spec)
  p_ab_pos <- annual_non_suppression * ab_pos_if_viremic + (1 - annual_non_suppression) * ab_pos_if_not
  pcr_confirm_per_year_triage <- p_ab_pos * ab_tests_per_year
  
  # cost per person-year under triage:
  cost_triage <- ab_tests_per_year * (cost_ab + cost_visit) + pcr_confirm_per_year_triage * (cost_pcr + cost_visit)
  
  # Effectiveness under triage: proportion of viremic detected earlier depends on testing interval and ab_sensitivity:
  # Assume mean delay to detection under periodic Ab testing = testing_interval/2 + any assay-specific delay (e.g., AB_rebound_delay). Here we assume detection interval = frequency_yrs/2.
  delay_triage_years <- frequency_yrs / 2
  
  # The earlier detection relative to PCR-only is: delay_pcronly_years - delay_triage_years
  # But only if triage delay < pcr delay, else no earlier detection.
  delta_delay <- pmax(0, delay_pcronly_years - delay_triage_years)
  
  # expected effect per person-year = (prob of being viremic in year)*(prob Ab detects them when tested during the period relative to PCR)*prob re-suppression * years gained
  # For simplicity: assume proportion of viremic detected by triage during the year = ab_sensitivity * (1 - delay_triage_years) ??? (too complex).
  # We'll approximate that testing at interval f will catch a fraction sens_frac = frequency_yrs/(frequency_yrs + delay_pcronly_years) * ab_sens
  # Simpler approach: suppose fraction of viremic caught earlier = ab_sens (since those who are viremic in the year are testable) * indicator(delta_delay>0)
  frac_caught_earlier <- ifelse(delta_delay > 0, ab_sens, 0)
  
  # Adjust probability of re-suppression depending on resistance/behaviour:
  # Weighted expected re-suppression probability among those detected earlier:
  p_resupp_weighted_earlier <- prop_res * p_resup_res * eff_mult_res +
    (1 - prop_res) * p_resup_beh * eff_mult_beh
  
  eff_triage <- annual_non_suppression * frac_caught_earlier * p_resupp_weighted_earlier * delta_delay
  
  # return per-person-year costs and effects
  tibble(
    ab_tests_per_year = ab_tests_per_year,
    pcr_tests_per_year_pcronly = pcr_tests_per_year_pcronly,
    pcr_tests_per_year_triage = pcr_confirm_per_year_triage,
    cost_pcronly = cost_pcronly,
    cost_triage = cost_triage,
    eff_pcronly = eff_pcronly,
    eff_triage = eff_triage,
    delta_cost = cost_triage - cost_pcronly,
    delta_eff = eff_triage - eff_pcronly,
    ICER = ifelse(delta_eff == 0, NA, delta_cost / delta_eff)
  )
}

##Specifying the frequencies that should be outputed
freq_scenarios <- tibble(
  frequency_yrs = c(1, 0.5, 0.25, 0.125),
  freq_label = c("Annual", "Biannual", "Quarterly", "6-weekly")
)

# Test function with some frequencies
calc_per_person(1, params)
calc_per_person(0.5, params)
calc_per_person(0.0833, params) # ~monthly

# ---------------------------------
# 3) Sweep frequencies & non-suppression prevalence
# ---------------------------------
freqs_yrs <- c(1, 0.5, 0.25, 0.125, 1/12)  # annual, biannual, quarterly, 6-week ~0.125yrs, monthly
prevalences <- c(0.02, 0.05, 0.10, 0.20)   # sensitivity to non-suppression rate

grid <- expand.grid(frequency_yrs = freqs_yrs, annual_non_suppression = prevalences) %>%
  arrange(annual_non_suppression, frequency_yrs) %>%
  as_tibble()

results <- grid %>%
  rowwise() %>%
  mutate(out = list(calc_per_person(frequency_yrs = frequency_yrs,
                                    params_list = params,
                                    annual_non_suppression = annual_non_suppression))) %>%
  unnest(out)

# Add labels and compute whether triage cheaper than PCR-only (Q2)
results <- results %>%
  mutate(
    freq_label = case_when(
      frequency_yrs == 1 ~ "annual",
      frequency_yrs == 0.5 ~ "biannual",
      frequency_yrs == 0.25 ~ "quarterly",
      frequency_yrs == 0.125 ~ "6-week",
      abs(frequency_yrs - 1/12) < 1e-6 ~ "monthly",
      TRUE ~ paste0(round(1/frequency_yrs,1),"x/yr")
    ),
    freq_label = str_trim(freq_label),
    freq_label = factor(
      freq_label,
      levels = c("6-weekly", "Quarterly", "Biannual", "Annual")
    ),
    triage_cheaper = cost_triage < cost_pcronly,
    triage_more_effective = eff_triage > eff_pcronly
  )

# Print results table
print(results %>% 
        dplyr::select(annual_non_suppression, freq_label, cost_pcronly, cost_triage, delta_cost, eff_pcronly, eff_triage, delta_eff, ICER, triage_cheaper, triage_more_effective))

# Answer to Q2: at what frequency cost of testing becomes lower than PCR-only?
# For each prevalence, find most frequent (largest ab_tests_per_year) that still has triage_cheaper == TRUE
q2_summary <- results %>%
  group_by(annual_non_suppression) %>%
  filter(triage_cheaper) %>%
  summarize(min_freq_yrs = max(frequency_yrs),   # largest interval (i.e., least frequent) where cheaper; adjust depending on interpretation
            freq_label = freq_label[which.max(frequency_yrs)],
            .groups = "drop")

print(q2_summary)

# ---------------------------------
# 4) Cost-effectiveness across WTP thresholds (Q3 + Q4)
# ---------------------------------
# For each frequency and prevalence compute net monetary benefit (NMB) = WTP*delta_eff - delta_cost
ce_results <- results %>%
  rowwise() %>%
  mutate(
    NMB_200 = WTPs[1] * delta_eff - delta_cost,
    NMB_500 = WTPs[2] * delta_eff - delta_cost,
    NMB_800 = WTPs[3] * delta_eff - delta_cost
  ) %>%
  ungroup()

# For each WTP and prevalence, find frequencies with positive NMB (i.e., triage cost-effective)
ce_summary <- ce_results %>%
  pivot_longer(cols = starts_with("NMB_"), names_to = "WTP_label", values_to = "NMB") %>%
  mutate(WTP = case_when(WTP_label == "NMB_200" ~ 200,
                         WTP_label == "NMB_500" ~ 500,
                         WTP_label == "NMB_800" ~ 800)) %>%
  group_by(annual_non_suppression, WTP) %>%
  summarize(
    best_freq = freq_label[which.max(NMB)],
    best_NMB = max(NMB),
    any_CE = any(NMB > 0),
    .groups = "drop"
  )

print(ce_summary)

# ---------------------------------
# 5) Tornado (univariate) sensitivity analysis
# ---------------------------------
# Define base and ranges for key inputs (examples)
tornado_inputs <- tibble(
  param = c("ab_sensitivity","ab_specificity","cost_ab","cost_pcr","mean_gain","prop_res","p_resup_res","p_resup_beh","annual_non_suppression"),
  base  = c(params$ab_sensitivity, params$ab_specificity, params$cost_ab, params$cost_pcr, params$mean_effect_years_gain_on_detection, params$prop_reistance, params$p_resupp_after_resistance_detect, params$p_resupp_after_behav_detect, params$annual_non_suppression),
  low   = c(0.3, 0.85, 1.0, 20, 0.05, 0.1, 0.05, 0.2, 0.01),
  high  = c(0.9, 0.99, 7.0, 100, 0.5, 0.7, 0.5, 0.95, 0.30)
)

# pick a scenario to test (e.g., biannual)
# base_freq <- 0.5
tornado_all <- freq_scenarios %>%
  rowwise() %>%
  mutate(
    tornado = list({
      base_freq <- frequency_yrs
      
      tornado_out <- tornado_inputs %>%
        rowwise() %>%
        mutate(
          out_low = list({
            pmod <- params
            pmod[[param]] <- low
            calc_per_person(base_freq, pmod)
          }),
          out_high = list({
            pmod <- params
            pmod[[param]] <- high
            calc_per_person(base_freq, pmod)
          })
        ) %>%
        unnest(cols = c(out_low, out_high), names_sep = "_") %>%
        mutate(
          NMB_low  = 500 * (out_low_eff_triage  - out_low_eff_pcronly) -
            (out_low_cost_triage     - out_low_cost_pcronly),
          NMB_high = 500 * (out_high_eff_triage - out_high_eff_pcronly) -
            (out_high_cost_triage    - out_high_cost_pcronly),
          NMB_base = {
            base <- calc_per_person(base_freq, params)
            500 * base$delta_eff - base$delta_cost
          },
          min_change = pmin(NMB_low - NMB_base, NMB_high - NMB_base),
          max_change = pmax(NMB_low - NMB_base, NMB_high - NMB_base)#,
          # freq_label = freq_label
        )
      
      tornado_out
    })
  ) %>%
  unnest(tornado)

tornado_plot_all <- ggplot(tornado_all) +
  geom_linerange(
    aes(x = param, ymin = min_change, ymax = max_change),
    size = 5
  ) +
  coord_flip() +
  # facet_wrap(~ freq_label, scales = "free_x") +
  facet_wrap(
    ~ forcats::fct_relevel(
      freq_label,
      "6-weekly",
      "Quarterly",
      "Biannual",
      "Annual"
    ),
    scales = "free_x"
  ) +
  labs(
    x = "",
    y = "Change in NMB (USD, WTP = 500)",
    title = "One-way sensitivity analysis across monitoring frequencies"
  ) +
  theme_minimal()

print(tornado_plot_all)

# ---------------------------------
# 6) Probabilistic Sensitivity Analysis (PSA)
# ---------------------------------
# Define distributions for uncertain parameters. Use Monte Carlo draws.
n_sims <- 5000

# helper to draw
draws <- tibble(
  ab_sensitivity = rbeta(n_sims, 60, 40),       # mean ~0.6 (alpha/beta tuned)
  ab_specificity = rbeta(n_sims, 95, 5),        # mean ~0.95
  cost_ab = rgamma(n_sims, shape = 4, scale = 0.6), # mean ~2.4
  cost_pcr = rgamma(n_sims, shape = 40, scale = 1.0), # mean ~40
  mean_gain = rlnorm(n_sims, meanlog = log(params$mean_effect_years_gain_on_detection), sdlog = 0.6),
  prop_res = rbeta(n_sims, 40, 60),             # mean ~0.4
  p_resup_res = rbeta(n_sims, 20, 80),          # mean ~0.2
  p_resup_beh = rbeta(n_sims, 60, 40),          # mean ~0.6
  annual_non_suppression = rbeta(n_sims, 10, 90) # mean ~0.1
)

# Run PSA for two policies at a chosen frequency (e.g., biannual)
# psa_freq <- 0.5
psa_all <- freq_scenarios %>%
  rowwise() %>%
  mutate(
    psa = list({
      psa_freq <- frequency_yrs
      
      draws %>%
        mutate(sim = row_number()) %>%
        rowwise() %>%
        mutate(
          out = list(calc_per_person(
            frequency_yrs = psa_freq,
            params_list = params,
            annual_non_suppression = annual_non_suppression,
            ab_sens = ab_sensitivity,
            ab_spec = ab_specificity,
            cost_ab = cost_ab,
            cost_pcr = cost_pcr,
            mean_gain = mean_gain,
            prop_res = prop_res,
            p_resup_res = p_resup_res,
            p_resup_beh = p_resup_beh
          ))
        ) %>%
        unnest(out) %>%
        transmute(
          delta_cost,
          delta_eff#,
          # freq_label = freq_label
        )
    })
  ) %>%
  unnest(psa)

ce_plane_all <- ggplot(psa_all, aes(x = delta_eff, y = delta_cost)) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # facet_wrap(~ freq_label) 
  facet_wrap(
    ~ forcats::fct_relevel(
      freq_label,
      "6-weekly",
      "Quarterly",
      "Biannual",
      "Annual"
    ),
    scales = "free_x"
  ) +
  labs(
    x = "Δ Effect (years of viral suppression)",
    y = "Δ Cost (USD)",
    title = "Cost-effectiveness planes by monitoring frequency"
  ) +
  theme_minimal()

print(ce_plane_all)

# CEAC: proportion of sims where NMB > 0 across WTPs
ceac_all <- psa_all %>%
  crossing(WTP = seq(0, 1000, by = 20)) %>%
  mutate(NMB = WTP * delta_eff - delta_cost) %>%
  group_by(freq_label, WTP) %>%
  summarize(
    pr_CE = mean(NMB > 0),
    .groups = "drop"
  )

ceac_plot_all <- ggplot(ceac_all, aes(x = WTP, y = pr_CE)) +
  geom_line(size = 1) +
  # facet_wrap(~ freq_label)
  facet_wrap(
    ~ forcats::fct_relevel(
      freq_label,
      "6-weekly",
      "Quarterly",
      "Biannual",
      "Annual"
    ),
    scales = "free_x"
  )+
  labs(
    x = "Willingness-to-pay (USD per year of suppression)",
    y = "Probability cost-effective",
    title = "CEACs by monitoring frequency"
  ) +
  theme_minimal()

print(ceac_plot_all)

# ---------------------------------
# 7) Outputs: Figures and summary tables
# ---------------------------------
# Save main results tables
write.csv(results, "CE_AB_analysis/ab_triage_results_sweep.csv", row.names = FALSE)
write.csv(ce_summary, "CE_AB_analysis/ab_triage_ce_summary.csv", row.names = FALSE)

# Save example plots
ggsave("CE_AB_analysis/tornado_all_frequencies.png", tornado_plot_all, width = 12, height = 8)
ggsave("CE_AB_analysis/ce_plane_all_frequencies.png", ce_plane_all, width = 12, height = 8)
ggsave("CE_AB_analysis/ceac_all_frequencies.png", ceac_plot_all, width = 12, height = 8)


message("Completed. Results saved to CSV and plots saved. Adjust parameters in 'params' and rerun as needed.")
