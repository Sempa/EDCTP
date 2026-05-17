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
library(patchwork)
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
  cost_visit_ab = 0,               # Cost of AB clinic visit if we use POC testing 
  
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
  efficacy_multiplier_behav = 1.0,          # detection more effective for behaviour
  # Mortality among suppressed and none suppressed
  mort_rate_suppressed = 0.015,
  mort_rate_unsuppressed = 0.06,
  # DALY parameters
  life_expectancy_remaining = 30,   # average remaining life-years at cohort age
  dw_suppressed = 0.05,             # disability weight (stable HIV on ART)
  dw_unsuppressed = 0.20            # disability weight (unsuppressed HIV)
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
                            cost_visit_ab = params_list$cost_clinic_visit,
                            cost_visit_pcr = params_list$cost_clinic_visit,
                            mean_gain = params_list$mean_effect_years_gain_on_detection,
                            prop_res = params_list$prop_reistance,
                            p_resup_res = params_list$p_resupp_after_resistance_detect,
                            p_resup_beh = params_list$p_resupp_after_behav_detect,
                            eff_mult_res = params_list$efficacy_multiplier_resistance,
                            eff_mult_beh = params_list$efficacy_multiplier_behav,
                            mort_rate_suppressed = params_list$mort_rate_suppressed,
                            mort_rate_unsuppressed = params_list$mort_rate_unsuppressed,
                            life_expectancy_remaining = params_list$life_expectancy_remaining,
                            dw_suppressed = params_list$dw_suppressed,
                            dw_unsuppressed = params_list$dw_unsuppressed
) {
  
  # ---- Mortality probabilities (annual) ----
  p_death_suppressed   <- 1 - exp(-mort_rate_suppressed)
  p_death_unsuppressed <- 1 - exp(-mort_rate_unsuppressed)
  
  # ---- Test frequency ----
  ab_tests_per_year <- ifelse(is.infinite(frequency_yrs) || frequency_yrs <= 0, 0, 1 / frequency_yrs)
  
  # ---- PCR-only strategy ----
  pcr_tests_per_year_pcronly <- 1  # routine annual PCR
  
  cost_pcronly <- pcr_tests_per_year_pcronly *
    (cost_pcr + cost_visit_pcr)
  
  # PCR-only: expected suppression-years regained (your existing shorthand)
  delay_pcronly_years <- 0.5 * (1 / pcr_tests_per_year_pcronly)
  
  p_resupp_weighted <- prop_res * p_resup_res + (1 - prop_res) * p_resup_beh
  eff_pcronly <- annual_non_suppression * p_resupp_weighted * mean_gain  # suppression-years regained per person-year
  
  # Convert eff_pcronly into an "effective" unsuppressed fraction (bounded to [0,1])
  unsuppressed_pcronly <- annual_non_suppression - eff_pcronly
  unsuppressed_pcronly <- pmin(1, pmax(0, unsuppressed_pcronly))
  
  # Mortality and DALY components for PCR-only (CONSISTENT with triage arm)
  mort_pcronly <- unsuppressed_pcronly * p_death_unsuppressed +
    (1 - unsuppressed_pcronly) * p_death_suppressed
  
  YLL_pcronly <- mort_pcronly * life_expectancy_remaining
  YLD_pcronly <- (1 - unsuppressed_pcronly) * dw_suppressed +
    unsuppressed_pcronly * dw_unsuppressed
  
  DALY_pcronly <- YLL_pcronly + YLD_pcronly
  
  
  # ---- AB-triage strategy ----
  # Ab test positive -> confirm PCR
  ab_pos_if_viremic <- ab_sens
  ab_pos_if_not     <- 1 - ab_spec
  
  p_ab_pos <- annual_non_suppression * ab_pos_if_viremic +
    (1 - annual_non_suppression) * ab_pos_if_not
  
  pcr_confirm_per_year_triage <- p_ab_pos * ab_tests_per_year
  
  cost_triage <- 
    ab_tests_per_year * (cost_ab + cost_visit_ab) +
    pcr_confirm_per_year_triage * (cost_pcr + cost_visit_pcr)
  
  # Detection timing
  delay_triage_years <- frequency_yrs / 2
  delta_delay <- pmax(0, delay_pcronly_years - delay_triage_years)
  
  frac_caught_earlier <- ifelse(delta_delay > 0, ab_sens, 0)
  
  p_resupp_weighted_earlier <- prop_res * p_resup_res * eff_mult_res +
    (1 - prop_res) * p_resup_beh * eff_mult_beh
  
  eff_triage <- annual_non_suppression * frac_caught_earlier * p_resupp_weighted_earlier * delta_delay
  
  # Convert eff_triage into an "effective" unsuppressed fraction (bounded to [0,1])
  unsuppressed_triage <- annual_non_suppression - eff_triage
  unsuppressed_triage <- pmin(1, pmax(0, unsuppressed_triage))
  
  mort_triage <- unsuppressed_triage * p_death_unsuppressed +
    (1 - unsuppressed_triage) * p_death_suppressed
  
  YLL_triage <- mort_triage * life_expectancy_remaining
  YLD_triage <- (1 - unsuppressed_triage) * dw_suppressed +
    unsuppressed_triage * dw_unsuppressed
  
  DALY_triage <- YLL_triage + YLD_triage
  
  
  # ---- Incremental outcomes ----
  delta_cost <- cost_triage - cost_pcronly
  
  # Suppression-years gained: triage vs PCR-only (as originally intended)
  delta_eff <- eff_triage - eff_pcronly
  
  # Life-years gained: convert deaths averted into LY using remaining life expectancy
  delta_deaths <- mort_pcronly - mort_triage
  delta_LY <- delta_deaths * life_expectancy_remaining
  
  # DALYs averted
  delta_DALY <- DALY_pcronly - DALY_triage
  
  
  # ---- ICERs ----
  ICER_suppression <- ifelse(delta_eff == 0, NA, delta_cost / delta_eff)
  ICER_LY          <- ifelse(delta_LY  == 0, NA, delta_cost / delta_LY)
  ICER_DALY        <- ifelse(delta_DALY== 0, NA, delta_cost / delta_DALY)
  
  tibble(
    ab_tests_per_year = ab_tests_per_year,
    pcr_tests_per_year_pcronly = pcr_tests_per_year_pcronly,
    pcr_tests_per_year_triage = pcr_confirm_per_year_triage,
    
    cost_pcronly = cost_pcronly,
    cost_triage  = cost_triage,
    
    eff_pcronly = eff_pcronly,
    eff_triage  = eff_triage,
    
    unsuppressed_pcronly = unsuppressed_pcronly,
    unsuppressed_triage  = unsuppressed_triage,
    
    mort_pcronly = mort_pcronly,
    mort_triage  = mort_triage,
    
    delta_cost = delta_cost,
    delta_eff  = delta_eff,
    delta_LY   = delta_LY,
    delta_DALY = delta_DALY,
    
    ICER_suppression = ICER_suppression,
    ICER_LY          = ICER_LY,
    ICER_DALY        = ICER_DALY,
    DALY_pcronly     = DALY_pcronly,
    DALY_triage      = DALY_triage
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

scenario_grid <- expand.grid(frequency_yrs = freqs_yrs, annual_non_suppression = prevalences) %>%
  arrange(annual_non_suppression, frequency_yrs) %>%
  as_tibble()

results <- scenario_grid %>%
  rowwise() %>%
  mutate(out = list(calc_per_person(frequency_yrs = frequency_yrs,
                                    params_list = params,
                                    annual_non_suppression = annual_non_suppression))) %>%
  unnest(out)

# Add labels and compute whether triage cheaper than PCR-only (Q2)
results <- results %>%
  mutate(
    freq_label = case_when(
      frequency_yrs == 1 ~ "Annual",
      frequency_yrs == 0.5 ~ "Biannual",
      frequency_yrs == 0.25 ~ "Quarterly",
      frequency_yrs == 0.125 ~ "6-weekly",
      abs(frequency_yrs - 1/12) < 1e-6 ~ "Monthly",
      TRUE ~ paste0(round(1/frequency_yrs,1),"x/yr")
    ),
    freq_label = factor(
      freq_label,
      levels = c("6-weekly", "Quarterly", "Biannual", "Annual", "Monthly")
    ),
    triage_cheaper = cost_triage < cost_pcronly,
    triage_more_effective = eff_triage > eff_pcronly
  )

# Print results table
print(results %>% 
        dplyr::select(annual_non_suppression, freq_label, cost_pcronly, 
                      cost_triage, delta_cost, eff_pcronly, eff_triage, delta_eff, 
                      ICER_suppression, ICER_LY, ICER_DALY, triage_cheaper, triage_more_effective))

# Answer to Q2: at what frequency cost of testing becomes lower than PCR-only?
# For each prevalence, find most frequent (largest ab_tests_per_year) that still has triage_cheaper == TRUE
q2_summary <- results %>%
  group_by(annual_non_suppression) %>%
  filter(triage_cheaper) %>%
  summarize(min_freq_yrs = max(frequency_yrs),   # largest interval (i.e., least frequent) where cheaper; adjust depending on interpretation
            freq_label = freq_label[which.max(frequency_yrs)],
            .groups = "drop")

print(q2_summary)

### Assuming we have POC AB testing and remove clinic visits
results_poc <- scenario_grid %>%
  
  rowwise() %>%
  
  mutate(
    out = list(
      calc_per_person(
        frequency_yrs = frequency_yrs,
        params_list = params,
        annual_non_suppression = annual_non_suppression,
        
        # Remove visit cost for Ab testing
        cost_visit_ab = 0
      )
    )
  ) %>%
  
  unnest(out) %>%
  
  mutate(
    
    # -----------------------------------------
    # Frequency labels
    # -----------------------------------------
    freq_label = case_when(
      frequency_yrs == 1     ~ "Annual",
      frequency_yrs == 0.5   ~ "Biannual",
      frequency_yrs == 0.25  ~ "Quarterly",
      frequency_yrs == 0.125 ~ "6-weekly",
      abs(frequency_yrs - 1/12) < 1e-6 ~ "Monthly",
      TRUE ~ paste0(round(1 / frequency_yrs, 1), "x/yr")
    ),
    
    freq_label = factor(
      freq_label,
      levels = c(
        "6-weekly",
        "Quarterly",
        "Biannual",
        "Annual",
        "Monthly"
      )
    ),
    
    # -----------------------------------------
    # Dominance indicators
    # -----------------------------------------
    triage_cheaper =
      cost_triage < cost_pcronly,
    
    triage_more_effective =
      eff_triage > eff_pcronly,
    
    # -----------------------------------------
    # Scenario label
    # -----------------------------------------
    scenario = "Point-of-care"
    
  )

results_facility <- results %>%
  mutate(
    scenario = "Facility-based"
  )

results_compare <- bind_rows(
  results_facility,
  results_poc
)

results_compare <- results_compare %>%
  mutate(
    NMB_500 = 500 * delta_DALY - delta_cost,
    prevalence_label = paste0(
      annual_non_suppression * 100,
      "% prevalence"
    )
  )

p_poc_nmb <- ggplot(
  results_compare,
  aes(
    x = freq_label,
    y = NMB_500,
    group = scenario,
    color = scenario
  )
) +
  
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  
  facet_wrap(~ prevalence_label) +
  
  scale_color_manual(
    values = c(
      "Facility-based" = "#d95f02",
      "Point-of-care" = "#1b9e77"
    )
  ) +
  
  labs(
    x = "Monitoring interval",
    y = "Net monetary benefit (US$)",
    color = "",
    title = "Economic impact of point-of-care antibody monitoring"
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

print(p_poc_nmb)

# ---------------------------------
# 4) Cost-effectiveness across WTP thresholds (Q3 + Q4)
# ---------------------------------
# For each frequency and prevalence compute net monetary benefit (NMB) = WTP*delta_eff - delta_cost
ce_results <- results %>%
  rowwise() %>%
  mutate(
    NMB_200 = WTPs[1] * delta_DALY - delta_cost,
    NMB_500 = WTPs[2] * delta_DALY - delta_cost,
    NMB_800 = WTPs[3] * delta_DALY - delta_cost
  ) %>%
  ungroup()

# For each WTP and prevalence, find frequencies with positive NMB (i.e., triage cost-effective)
ce_summary <- ce_results %>%
  group_by(annual_non_suppression, freq_label) %>%
  summarise(
    Cost_PCR_Only_USD = first(cost_pcronly),
    Cost_Ab_Triage_USD = first(cost_triage),
    Incremental_Cost_USD = first(delta_cost),
    
    Suppression_Years_Gained = first(delta_eff),
    Life_Years_Gained = first(delta_LY),
    DALYs_Averted = first(delta_DALY),
    
    ICER_Suppression = first(ICER_suppression),
    ICER_LY = first(ICER_LY),
    ICER_DALY = first(ICER_DALY),
    
    .groups = "drop"
  ) %>%
  arrange(annual_non_suppression, freq_label)

publication_table <- ce_summary %>%
  mutate(
    Cost_PCR_Only_USD = round(Cost_PCR_Only_USD, 2),
    Cost_Ab_Triage_USD = round(Cost_Ab_Triage_USD, 2),
    Incremental_Cost_USD = round(Incremental_Cost_USD, 2),
    
    Suppression_Years_Gained = round(Suppression_Years_Gained, 4),
    Life_Years_Gained = round(Life_Years_Gained, 4),
    DALYs_Averted = round(DALYs_Averted, 4),
    
    ICER_Suppression = round(ICER_Suppression, 2),
    ICER_LY = round(ICER_LY, 2),
    ICER_DALY = round(ICER_DALY, 2),
    Interpretation = case_when(
      Incremental_Cost_USD < 0 & DALYs_Averted > 0 ~ "Dominant",
      Incremental_Cost_USD > 0 & DALYs_Averted < 0 ~ "Dominated",
      ICER_DALY < 500 ~ "Cost-effective",
      TRUE ~ "Not cost-effective"
    )
  ) %>%
  rename(
    `Viral non-suppression prevalence` = annual_non_suppression,
    `Monitoring interval` = freq_label,
    `Annual cost: PCR-only (US$)` = Cost_PCR_Only_USD,
    `Annual cost: Ab-triage (US$)` = Cost_Ab_Triage_USD,
    `Incremental cost (US$)` = Incremental_Cost_USD,
    `Suppression-years gained` = Suppression_Years_Gained,
    `Life-years gained` = Life_Years_Gained,
    `DALYs averted` = DALYs_Averted,
    `ICER (US$/suppression-year)` = ICER_Suppression,
    `ICER (US$/life-year)` = ICER_LY,
    `ICER (US$/DALY averted)` = ICER_DALY,
    Interpretation = Interpretation
  )

print(publication_table)

# ============================================================
# TABLE: Absolute outcomes for PCR-only versus Ab-triage
# ============================================================

absolute_table <- results %>%
  
  mutate(
    
    # -----------------------------
    # Human-readable monitoring labels
    # -----------------------------
    `Monitoring interval` = freq_label,
    
    # -----------------------------
    # PCR-only strategy outcomes
    # -----------------------------
    `Annual cost: PCR-only (US$)` =
      round(cost_pcronly, 2),
    
    `Suppression-years: PCR-only` =
      round(eff_pcronly, 4),
    
    `Mean time unsuppressed: PCR-only (years)` =
      round(unsuppressed_pcronly, 4),
    
    `Mortality risk: PCR-only` =
      round(mort_pcronly, 4),
    
    `DALYs: PCR-only` =
      round(DALY_pcronly, 4),
    
    # -----------------------------
    # Antibody-triage strategy outcomes
    # -----------------------------
    `Annual cost: Ab-triage (US$)` =
      round(cost_triage, 2),
    
    `Suppression-years: Ab-triage` =
      round(eff_triage, 4),
    
    `Mean time unsuppressed: Ab-triage (years)` =
      round(unsuppressed_triage, 4),
    
    `Mortality risk: Ab-triage` =
      round(mort_triage, 4),
    
    `DALYs: Ab-triage` =
      round(DALY_triage, 4)
    
  ) %>%
  
  # -----------------------------
# Keep only variables for publication
# -----------------------------
dplyr::select(
  `Monitoring interval`,
  
  `Annual cost: PCR-only (US$)`,
  `Suppression-years: PCR-only`,
  `Mean time unsuppressed: PCR-only (years)`,
  `Mortality risk: PCR-only`,
  `DALYs: PCR-only`,
  
  `Annual cost: Ab-triage (US$)`,
  `Suppression-years: Ab-triage`,
  `Mean time unsuppressed: Ab-triage (years)`,
  `Mortality risk: Ab-triage`,
  `DALYs: Ab-triage`
)

# View table
print(absolute_table)
# ---------------------------------
# 5) Tornado (univariate) sensitivity analysis
# ---------------------------------
# Define base and ranges for key inputs (examples)
tornado_inputs <- tibble(
  param = c(
    "ab_sensitivity","ab_specificity",
    "cost_ab","cost_pcr",
    "mean_gain","prop_res",
    "p_resup_res","p_resup_beh",
    "annual_non_suppression",
    "mort_rate_suppressed",
    "mort_rate_unsuppressed"
  ),
  base  = c(
    params$ab_sensitivity,
    params$ab_specificity,
    params$cost_ab,
    params$cost_pcr,
    params$mean_effect_years_gain_on_detection,
    params$prop_reistance,
    params$p_resupp_after_resistance_detect,
    params$p_resupp_after_behav_detect,
    params$annual_non_suppression,
    params$mort_rate_suppressed,
    params$mort_rate_unsuppressed
  ),
  low   = c(
    0.3, 0.85,
    1.0, 20,
    0.05, 0.1,
    0.05, 0.2,
    0.01,
    0.005,   # suppressed mortality lower bound
    0.03     # unsuppressed lower bound
  ),
  high  = c(
    0.9, 0.99,
    7.0, 100,
    0.5, 0.7,
    0.5, 0.95,
    0.30,
    0.03,    # suppressed upper bound
    0.10     # unsuppressed upper bound
  )
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
          NMB_low  = 500 * out_low_delta_DALY  - out_low_delta_cost,
          NMB_high = 500 * out_high_delta_DALY - out_high_delta_cost,
          NMB_base = {
            base <- calc_per_person(base_freq, params)
            500 * base$delta_DALY - base$delta_cost
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
    y = "Change in NMB (USD, WTP = 500 per DALY averted)",
    title = "One-way sensitivity analysis across monitoring frequencies"
  ) +
  theme_minimal()

print(tornado_plot_all)

deterministic_plot_data <- results %>%
  mutate(
    prevalence_label = paste0(round(annual_non_suppression * 100), "%"),
    NMB_500 = 500 * delta_DALY - delta_cost
  ) %>%
  select(
    freq_label,
    prevalence_label,
    delta_cost,
    delta_DALY,
    ICER_DALY,
    NMB_500
  ) %>%
  pivot_longer(
    cols = c(delta_cost, delta_DALY, ICER_DALY, NMB_500),
    names_to = "Outcome",
    values_to = "Value"
  ) %>%
  mutate(
    Outcome = recode(
      Outcome,
      delta_cost = "Incremental cost (US$)",
      delta_DALY = "DALYs averted",
      ICER_DALY = "ICER (US$/DALY averted)",
      NMB_500 = "Net monetary benefit (US$)"
    )
  )


plot_cost <- deterministic_plot_data %>%
  filter(Outcome == "Incremental cost (US$)") %>%
  ggplot(aes(freq_label, Value, group = prevalence_label, color = prevalence_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(x = "Monitoring interval", y = "US$", color = "Prevalence") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

plot_daly <- deterministic_plot_data %>%
  filter(Outcome == "DALYs averted") %>%
  ggplot(aes(freq_label, Value, group = prevalence_label, color = prevalence_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(x = "Monitoring interval", y = "DALYs", color = "Prevalence") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

plot_icer <- deterministic_plot_data %>%
  filter(Outcome == "ICER (US$/DALY averted)") %>%
  ggplot(aes(freq_label, Value, group = prevalence_label, color = prevalence_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(x = "Monitoring interval", y = "US$/DALY", color = "Prevalence") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

plot_nmb <- deterministic_plot_data %>%
  filter(Outcome == "Net monetary benefit (US$)") %>%
  ggplot(aes(freq_label, Value, group = prevalence_label, color = prevalence_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(x = "Monitoring interval", y = "US$", color = "Prevalence") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

final_deterministic_figure <- (plot_cost + plot_daly) /
  (plot_icer + plot_nmb) +
  plot_annotation(
    tag_levels = "A"#,
    # title = "Deterministic cost-effectiveness of antibody triage by monitoring frequency"
  )

print(final_deterministic_figure)

# ---------------------------------
# 6) Probabilistic Sensitivity Analysis (PSA)
# ---------------------------------
# Define distributions for uncertain parameters. Use Monte Carlo draws.
n_sims <- 5000

# helper to draw
draws <- tibble(
  ab_sensitivity = rbeta(n_sims, 60, 40),
  ab_specificity = rbeta(n_sims, 95, 5),
  
  cost_ab = rgamma(n_sims, shape = 4, scale = 0.6),
  cost_pcr = rgamma(n_sims, shape = 40, scale = 1.0),
  
  mean_gain = rlnorm(
    n_sims,
    meanlog = log(params$mean_effect_years_gain_on_detection),
    sdlog = 0.6
  ),
  
  prop_res = rbeta(n_sims, 40, 60),
  p_resup_res = rbeta(n_sims, 20, 80),
  p_resup_beh = rbeta(n_sims, 60, 40),
  
  annual_non_suppression = rbeta(n_sims, 10, 90),
  mort_rate_suppressed = rbeta(n_sims, 2, 198),   # mean ~0.01
  mort_rate_unsuppressed = rbeta(n_sims, 6, 94)   # mean ~0.06
)

# Run PSA for two policies at a chosen frequency (e.g., biannual)
# psa_freq <- 0.5
WTP_threshold <- 500
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
            p_resup_beh = p_resup_beh,
            mort_rate_suppressed = mort_rate_suppressed,
            mort_rate_unsuppressed = mort_rate_unsuppressed
          ))
        ) %>%
        unnest(out) %>%
        transmute(#sim, 
          delta_cost, 
          delta_eff, 
          delta_LY, delta_DALY)
    })
  ) %>%
  unnest(psa) %>%
  mutate(
    delta_days = delta_DALY * 365,
    CE_status = ifelse(delta_cost <= WTP_threshold * delta_DALY,
                       "Cost-effective",
                       "Not cost-effective")
  ) %>%
  mutate(
    panel_label = case_when(
      freq_label == "6-weekly" ~ "A. Every 6 weeks",
      freq_label == "Quarterly" ~ "B. Every 3 months",
      freq_label == "Biannual" ~ "C. Every 6 months",
      freq_label == "Annual" ~ "D. Every 12 months"
    )
  )

ce_plane_all <- ggplot(
  psa_all,
  aes(x = delta_days, y = delta_cost, color = CE_status)
) +
  geom_point(alpha = 0.35, size = 1.2) +
  
  geom_abline(
    intercept = 0,
    slope = WTP_threshold / 365,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.4) +
  
  facet_wrap(
    ~ forcats::fct_relevel(
      freq_label,
      "6-weekly",
      "Quarterly",
      "Biannual",
      "Annual"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Cost-effective" = "#1b9e77",
      "Not cost-effective" = "#d95f02"
    )
  ) +
  
  labs(
    x = "Health gain (DALYs averted, days)",
    y = "Additional cost (US$)",
    color = "Simulation result"#,
    # title = "Probabilistic cost-effectiveness of antibody-guided monitoring by testing interval"
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  facet_wrap(~panel_label)

print(ce_plane_all)

# print(ce_plane_all)

# CEAC: proportion of sims where NMB > 0 across WTPs
ceac_all <- psa_all %>%
  crossing(WTP = seq(0, 1000, by = 20)) %>%
  mutate(NMB = WTP * delta_DALY - delta_cost) %>%
  group_by(freq_label, WTP) %>%
  summarize(
    pr_CE = mean(NMB > 0),
    .groups = "drop"
  ) %>%
  mutate(
    panel_label = case_when(
      freq_label == "6-weekly" ~ "A. Every 6 weeks",
      freq_label == "Quarterly" ~ "B. Every 3 months",
      freq_label == "Biannual" ~ "C. Every 6 months",
      freq_label == "Annual" ~ "D. Every 12 months"
    ),
    pr_CE_percent = pr_CE * 100
  )

ceac_plot_all <- ggplot(
  ceac_all,
  aes(x = WTP, y = pr_CE_percent)
) +
  geom_line(linewidth = 1.2) +
  
  facet_wrap(
    ~ factor(
      panel_label,
      levels = c(
        "A. Every 6 weeks",
        "B. Every 3 months",
        "C. Every 6 months",
        "D. Every 12 months"
      )
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    labels = function(x) paste0(x, "%")
  ) +
  
  labs(
    x = "Willingness-to-pay threshold (US$ per DALY averted)",
    y = "Probability the strategy is cost-effective"#,
    # title = "Probability that antibody-guided monitoring is cost-effective at different testing intervals"
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

print(ceac_plot_all)


# ---------------------------------
# 7) Outputs: Figures and summary tables
# ---------------------------------
# Save main results tables
write.csv(results, "CE_AB_analysis/ab_triage_results_sweep.csv", row.names = FALSE)
write.csv(publication_table, "CE_AB_analysis/ab_triage_publication_table.csv", row.names = FALSE)
write.csv(absolute_table, "CE_AB_analysis/absolute_table_publication_table.csv", row.names = FALSE)
# Save example plots
ggsave("CE_AB_analysis/tornado_all_frequencies.png", tornado_plot_all, width = 12, height = 8)
ggsave("CE_AB_analysis/deterministic_headline_figure.png", final_deterministic_figure, width = 12, height = 8, dpi = 600)
ggsave("CE_AB_analysis/ce_plane_all_frequencies.png", ce_plane_all, width = 12, height = 8)
ggsave("CE_AB_analysis/ceac_all_frequencies.png", ceac_plot_all, width = 12, height = 8)
ggsave("CE_AB_analysis/p_poc_nmb.png", p_poc_nmb, width = 10, height = 8)

message("Completed. Results saved to CSV and plots saved. Adjust parameters in 'params' and rerun as needed.")
