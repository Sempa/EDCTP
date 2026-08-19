# ==============================================================================
# Master Cost-Effectiveness Analysis Script: Ab-Triage vs PCR & Hybrid
# Integrated and Refined Architecture
# Requirements: tidyverse, data.table, ggplot2, cowplot, reshape2, scales, patchwork, gt, knitr
# ==============================================================================

library(tidyverse)
library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(patchwork)
library(gt)
library(knitr)

set.seed(12345)

# ------------------------------------------------------------------------------
# 1) Global Model Parameters & Thresholds
# ------------------------------------------------------------------------------
params <- list(
  # Population & epidemic
  N = 10000,                       # cohort size
  annual_non_suppression = 0.10,   # baseline annual proportion with VL non-suppression
  
  # Test characteristics
  ab_sensitivity = 0.60,           # per-test sensitivity of Ab test
  ab_specificity = 0.95,           # specificity of Ab test
  pcr_sensitivity = 1.0,           # assumed
  pcr_specificity = 1.0,           # assumed
  antibody_lag = 0.08,             # ~4 weeks lag between viremia and detectable Ab
  
  # Costs (USD)
  cost_ab = 2.0,                   # cost per Ab test
  cost_pcr = 40.0,                 # cost per PCR test
  cost_clinic_visit = 5.0,         # clinic visit cost per test round
  cost_visit_ab = 5.0,             # baseline facility cost (0 for POC)
  
  # Effectiveness & Resuppression parameters
  mean_gain_delay_reduction = 0.25, # baseline gain factor
  prop_resistance = 0.4,           # proportion of viremic with resistance
  p_resupp_after_resistance_detect = 0.2,
  p_resupp_after_behav_detect = 0.6,
  efficacy_multiplier_resistance = 0.6,
  efficacy_multiplier_behav = 1.0,
  
  # Mortality & DALYs
  mort_rate_suppressed = 0.015,
  mort_rate_unsuppressed = 0.06,
  life_expectancy_remaining = 30,
  dw_suppressed = 0.05,            # disability weight (suppressed)
  dw_unsuppressed = 0.20           # disability weight (unsuppressed)
)

WTPs <- c(200, 500, 800)

freq_scenarios <- tibble(
  frequency_yrs = c(1, 0.5, 0.25, 0.125),
  freq_label = c("Annual", "Biannual", "Quarterly", "6-weekly")
)

# Output Directory
dir.create("CE_AB_analysis", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2) Core Calculation Helper Function
# ------------------------------------------------------------------------------
calc_per_person <- function(frequency_yrs = 1,
                            strategy = "AB_triage", 
                            params_list = params,
                            annual_non_suppression = params_list$annual_non_suppression,
                            ab_sens = params_list$ab_sensitivity,
                            ab_spec = params_list$ab_specificity,
                            antibody_lag = params_list$antibody_lag,
                            cost_ab = params_list$cost_ab,
                            cost_pcr = params_list$cost_pcr,
                            cost_visit_ab = params_list$cost_visit_ab,
                            cost_visit_pcr = params_list$cost_clinic_visit,
                            prop_res = params_list$prop_resistance,
                            p_resup_res = params_list$p_resupp_after_resistance_detect,
                            p_resup_beh = params_list$p_resupp_after_behav_detect,
                            eff_mult_res = params_list$efficacy_multiplier_resistance,
                            eff_mult_beh = params_list$efficacy_multiplier_behav,
                            mort_rate_suppressed = params_list$mort_rate_suppressed,
                            mort_rate_unsuppressed = params_list$mort_rate_unsuppressed,
                            life_expectancy_remaining = params_list$life_expectancy_remaining,
                            dw_suppressed = params_list$dw_suppressed,
                            dw_unsuppressed = params_list$dw_unsuppressed,
                            wtp_threshold = 500) {
  
  # Rates & resuppression probability
  p_death_suppressed   <- 1 - exp(-mort_rate_suppressed)
  p_death_unsuppressed <- 1 - exp(-mort_rate_unsuppressed)
  
  p_resupp_weighted <- prop_res * p_resup_res * eff_mult_res + (1 - prop_res) * p_resup_beh * eff_mult_beh
  p_fail_resuppress <- 1 - p_resupp_weighted
  
  # PCR-only Baseline (Standard of Care)
  cost_pcronly <- 1 * (cost_pcr + cost_visit_pcr)
  delay_pcronly <- 0.5
  time_unsupp_pcronly <- annual_non_suppression * (delay_pcronly + (1 - delay_pcronly) * p_fail_resuppress)
  eff_pcronly <- annual_non_suppression - time_unsupp_pcronly
  
  ab_tests_per_year <- ifelse(is.infinite(frequency_yrs) || frequency_yrs <= 0, 0, 1 / frequency_yrs)
  interval_length <- frequency_yrs
  
  if (strategy == "PCR_only") {
    cost_strat <- cost_pcronly
    time_unsupp_strat <- time_unsupp_pcronly
    pcr_confirm_per_year <- 1
    ab_tests_count <- 0
    
  } else if (strategy %in% c("AB_triage", "Hybrid")) {
    post_lag_window <- pmax(0, 1.0 - antibody_lag)
    n_post_lag_tests <- pmax(0, post_lag_window / interval_length)
    cum_detection <- ifelse(n_post_lag_tests <= 0, 0, 1 - (1 - ab_sens)^n_post_lag_tests)
    delay_triage <- (interval_length / 2) + antibody_lag
    
    if (strategy == "AB_triage") {
      time_unsupp_strat <- annual_non_suppression * (
        cum_detection * (delay_triage + (1 - delay_triage) * p_fail_resuppress) +
          (1 - cum_detection) * (1.0)
      )
      
      p_ab_pos <- annual_non_suppression * ab_sens + (1 - annual_non_suppression) * (1 - ab_spec)
      ab_tests_count <- ab_tests_per_year
      pcr_confirm_per_year <- p_ab_pos * ab_tests_per_year
      
      cost_strat <- ab_tests_count * (cost_ab + cost_visit_ab) + 
        pcr_confirm_per_year * (cost_pcr + cost_visit_pcr)
      
    } else if (strategy == "Hybrid") {
      time_unsupp_strat <- annual_non_suppression * (
        cum_detection * (delay_triage + (1 - delay_triage) * p_fail_resuppress) +
          (1 - cum_detection) * (0.5 + 0.5 * p_fail_resuppress)
      )
      
      p_ab_pos <- annual_non_suppression * ab_sens + (1 - annual_non_suppression) * (1 - ab_spec)
      ab_tests_count <- ab_tests_per_year
      pcr_confirm_per_year <- 1 + (p_ab_pos * ab_tests_count)
      
      cost_strat <- ab_tests_count * (cost_ab + cost_visit_ab) + 
        pcr_confirm_per_year * (cost_pcr + cost_visit_pcr)
    }
  }
  
  time_unsupp_strat <- pmin(1, pmax(0, time_unsupp_strat))
  eff_triage <- annual_non_suppression - time_unsupp_strat
  
  # Mortality and DALYs
  mort_strat <- time_unsupp_strat * p_death_unsuppressed + (1 - time_unsupp_strat) * p_death_suppressed
  YLL_strat  <- mort_strat * life_expectancy_remaining
  YLD_strat  <- (1 - time_unsupp_strat) * dw_suppressed + time_unsupp_strat * dw_unsuppressed
  DALY_strat <- YLL_strat + YLD_strat
  
  mort_pcronly <- time_unsupp_pcronly * p_death_unsuppressed + (1 - time_unsupp_pcronly) * p_death_suppressed
  DALY_pcronly <- (mort_pcronly * life_expectancy_remaining) + 
    ((1 - time_unsupp_pcronly) * dw_suppressed + time_unsupp_pcronly * dw_unsuppressed)
  
  # Incrementals
  delta_cost <- cost_strat - cost_pcronly
  delta_eff  <- eff_triage - eff_pcronly
  delta_LY   <- (mort_pcronly - mort_strat) * life_expectancy_remaining
  delta_DALY <- DALY_pcronly - DALY_strat
  nmb        <- wtp_threshold * delta_DALY - delta_cost
  
  # ICERs
  ICER_suppression <- ifelse(delta_eff == 0, NA, delta_cost / delta_eff)
  ICER_LY          <- ifelse(delta_LY == 0, NA, delta_cost / delta_LY)
  ICER_DALY        <- ifelse(delta_DALY == 0, NA, delta_cost / delta_DALY)
  
  verdict <- case_when(
    delta_cost < 0 & delta_DALY > 0 ~ "Dominant",
    nmb > 0 & delta_DALY > 0        ~ "Cost-effective",
    delta_DALY > 0                  ~ "Health-improving (Not CE)",
    TRUE                            ~ "Health-harming"
  )
  
  tibble(
    strategy = strategy,
    ab_tests_per_year = ab_tests_count,
    pcr_tests_per_year_pcronly = 1,
    pcr_tests_per_year_triage = pcr_confirm_per_year,
    cost = cost_strat,
    cost_pcronly = cost_pcronly,
    cost_triage = cost_strat,
    eff_pcronly = eff_pcronly,
    eff_triage = eff_triage,
    unsuppressed_pcronly = time_unsupp_pcronly,
    unsuppressed_triage = time_unsupp_strat,
    mort_pcronly = mort_pcronly,
    mort_triage = mort_strat,
    DALY_pcronly = DALY_pcronly,
    DALY_triage = DALY_strat,
    time_unsuppressed = time_unsupp_strat,
    DALY = DALY_strat,
    delta_cost = delta_cost,
    delta_eff = delta_eff,
    delta_LY = delta_LY,
    delta_DALY = delta_DALY,
    ICER_suppression = ICER_suppression,
    ICER_LY = ICER_LY,
    ICER_DALY = ICER_DALY,
    nmb = nmb,
    verdict = verdict
  )
}

# ------------------------------------------------------------------------------
# 3) Decision Map: Sensitivity vs Frequency
# ------------------------------------------------------------------------------
grid_decision <- expand.grid(
  ab_sens = seq(0.40, 0.99, length.out = 100),
  freq_num = seq(1, 12, length.out = 100),
  strategy = "AB_triage"
) %>% 
  as_tibble() %>% 
  mutate(frequency_yrs = 1 / freq_num)

results_decision <- grid_decision %>%
  rowwise() %>%
  mutate(out = list(calc_per_person(
    frequency_yrs = frequency_yrs, 
    strategy = strategy, 
    ab_sens = ab_sens
  ))) %>%
  select(-strategy) %>% # Drop the duplicate column before unnesting
  unnest(out) %>%
  mutate(verdict = factor(
    verdict, 
    levels = c("Dominant", "Cost-effective", "Health-improving (Not CE)", "Health-harming")
  ))

fig_decision_map <- ggplot(results_decision, aes(x = ab_sens, y = freq_num, fill = verdict)) +
  geom_tile() +
  scale_fill_manual(
    values = c("Dominant" = "#2b83ba", "Cost-effective" = "#abdda4", "Health-improving (Not CE)" = "#fdae61", "Health-harming" = "#d7191c"),
    drop = FALSE
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_y_continuous(
    breaks = c(1, 2, 4, 8, 12),
    labels = c("1 (Annual)", "2 (Biannual)", "4 (Quarterly)", "8 (6-Weekly)", "12 (Monthly)"),
    expand = c(0, 0)
  ) +
  labs(
    x = "Antibody Test Sensitivity (Per Test)",
    y = "Testing Frequency (Tests per Year)",
    fill = "Policy Verdict",
    title = "Operational Requirements for Antibody Triage Success",
    subtitle = "Decision boundaries across antibody sensitivity and monitoring frequency (WTP = $500/DALY)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid = element_blank(), plot.title = element_text(face = "bold"))

print(fig_decision_map)
ggsave("CE_AB_analysis/decision_map_sensitivity_vs_frequency.png", fig_decision_map, width = 10, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 4) Strategy Thresholds Table & Headline Figure
# ------------------------------------------------------------------------------
calculate_strategy_thresholds <- function(params_list = params, wtp = 500) {
  sweep_grid <- expand.grid(
    ab_sens = seq(0.40, 0.999, by = 0.001),
    frequency_yrs = c(0.5, 0.25, 0.125, 1/12),
    strategy = c("AB_triage", "Hybrid")
  ) %>% as_tibble()
  
  full_results <- sweep_grid %>%
    rowwise() %>%
    mutate(out = list(calc_per_person(
      frequency_yrs = frequency_yrs, 
      strategy = strategy, 
      ab_sens = ab_sens, 
      params_list = params_list, 
      wtp_threshold = wtp
    ))) %>%
    select(-strategy) %>% # Drop duplicate strategy column before unnesting
    unnest(out) %>%
    mutate(
      interval_label = case_when(
        frequency_yrs == 0.5 ~ "Biannual",
        frequency_yrs == 0.25 ~ "Quarterly",
        frequency_yrs == 0.125 ~ "6-Weekly",
        abs(frequency_yrs - 1/12) < 1e-5 ~ "Monthly"
      ),
      interval_label = factor(interval_label, levels = c("Biannual", "Quarterly", "6-Weekly", "Monthly"))
    )
  
  full_results %>%
    group_by(strategy, interval_label) %>%
    summarise(
      min_sens_health_gain = ab_sens[which(delta_DALY > 0)[1]],
      min_sens_cost_effective = ab_sens[which(nmb > 0 & delta_DALY > 0)[1]],
      min_sens_dominant = ab_sens[which(delta_DALY > 0 & delta_cost < 0)[1]],
      .groups = "drop"
    ) %>%
    mutate(
      min_sens_health_gain = ifelse(is.na(min_sens_health_gain), "Not Viable", scales::percent(min_sens_health_gain, accuracy = 0.1)),
      min_sens_cost_effective = ifelse(is.na(min_sens_cost_effective), "Not Viable", scales::percent(min_sens_cost_effective, accuracy = 0.1)),
      min_sens_dominant = ifelse(is.na(min_sens_dominant), "Not Viable", scales::percent(min_sens_dominant, accuracy = 0.1))
    )
}

threshold_summary <- calculate_strategy_thresholds(params, wtp = 500)

gt_threshold_table <- threshold_summary %>%
  gt(groupname_col = "strategy") %>%
  tab_header(
    title = md("**Minimum Required Antibody Sensitivity Thresholds**"),
    subtitle = "Assessing operational viability across testing frequencies and monitoring strategies"
  ) %>%
  cols_label(
    interval_label = "Monitoring Frequency",
    min_sens_health_gain = "Health Improvement (ΔDALY > 0)",
    min_sens_cost_effective = "Cost-Effective (NMB > 0)",
    min_sens_dominant = "Dominant (ΔDALY > 0 & ΔCost < 0)"
  ) %>%
  cols_align(align = "center", columns = c(min_sens_health_gain, min_sens_cost_effective, min_sens_dominant)) %>%
  opt_row_striping()

print(gt_threshold_table)

sens_sweep <- expand.grid(
  ab_sens = seq(0.40, 0.99, by = 0.005),
  frequency_yrs = c(0.5, 0.25, 0.125, 1/12),
  strategy = c("AB_triage", "Hybrid")
) %>% as_tibble()

sens_results <- sens_sweep %>%
  rowwise() %>%
  mutate(out = list(calc_per_person(frequency_yrs = frequency_yrs, strategy = strategy, ab_sens = ab_sens))) %>%
  select(-strategy) %>%
  unnest(out) %>%
  mutate(
    freq_label = case_when(
      frequency_yrs == 0.5 ~ "Biannual",
      frequency_yrs == 0.25 ~ "Quarterly",
      frequency_yrs == 0.125 ~ "6-weekly",
      abs(frequency_yrs - 1/12) < 1e-6 ~ "Monthly"
    ),
    freq_label = factor(freq_label, levels = c("Monthly", "6-weekly", "Quarterly", "Biannual"))
  )

fig_headline_threshold <- ggplot(sens_results, aes(x = ab_sens, y = delta_DALY, color = freq_label, linetype = strategy)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Antibody Test Sensitivity (Per Test)",
    y = "DALYs Averted per Person-Year vs PCR-Only Baseline",
    color = "Testing Interval",
    linetype = "Strategy",
    title = "Sensitivity Thresholds Required for Health Improvements",
    subtitle = "Accounting for 4-week seroconversion lag and continuous rebound arrival"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

print(fig_headline_threshold)
ggsave("CE_AB_analysis/headline_sensitivity_thresholds.png", fig_headline_threshold, width = 10, height = 6, dpi = 300)
write.csv(sens_results, "CE_AB_analysis/ab_sensitivity_threshold_sweep.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 5) Deterministic Scenario Sweeps & Publication Tables
# ------------------------------------------------------------------------------
freqs_yrs <- c(1, 0.5, 0.25, 0.125, 1/12)
prevalences <- c(0.02, 0.05, 0.10, 0.20)

scenario_grid <- expand.grid(
  frequency_yrs = freqs_yrs, 
  annual_non_suppression = prevalences,
  strategy = "AB_triage"
) %>% arrange(annual_non_suppression, frequency_yrs) %>% as_tibble()

results <- scenario_grid %>%
  rowwise() %>%
  mutate(out = list(calc_per_person(
    frequency_yrs = frequency_yrs, 
    strategy = strategy, 
    params_list = params, 
    annual_non_suppression = annual_non_suppression
  ))) %>%
  select(-strategy) %>% # Drop duplicate strategy column before unnesting
  unnest(out) %>%
  mutate(
    freq_label = case_when(
      frequency_yrs == 1 ~ "Annual",
      frequency_yrs == 0.5 ~ "Biannual",
      frequency_yrs == 0.25 ~ "Quarterly",
      frequency_yrs == 0.125 ~ "6-weekly",
      abs(frequency_yrs - 1/12) < 1e-6 ~ "Monthly",
      TRUE ~ paste0(round(1/frequency_yrs, 1), "x/yr")
    ),
    freq_label = factor(freq_label, levels = c("6-weekly", "Quarterly", "Biannual", "Annual", "Monthly")),
    triage_cheaper = cost_triage < cost_pcronly,
    triage_more_effective = eff_triage > eff_pcronly
  )

write.csv(results, "CE_AB_analysis/ab_triage_results_sweep.csv", row.names = FALSE)

q2_summary <- results %>%
  group_by(annual_non_suppression) %>%
  filter(triage_cheaper) %>%
  summarize(min_freq_yrs = max(frequency_yrs), freq_label = freq_label[which.max(frequency_yrs)], .groups = "drop")
print(q2_summary)

# Publication Tables
publication_table <- results %>%
  mutate(
    Cost_PCR_Only_USD = round(cost_pcronly, 2),
    Cost_Ab_Triage_USD = round(cost_triage, 2),
    Incremental_Cost_USD = round(delta_cost, 2),
    Suppression_Years_Gained = round(delta_eff, 4),
    Life_Years_Gained = round(delta_LY, 4),
    DALYs_Averted = round(delta_DALY, 4),
    ICER_Suppression = round(ICER_suppression, 2),
    ICER_LY = round(ICER_LY, 2),
    ICER_DALY = round(ICER_DALY, 2),
    Interpretation = case_when(
      Incremental_Cost_USD < 0 & DALYs_Averted > 0 ~ "Dominant",
      Incremental_Cost_USD > 0 & DALYs_Averted < 0 ~ "Dominated",
      ICER_DALY < 500 ~ "Cost-effective",
      TRUE ~ "Not cost-effective"
    )
  ) %>%
  dplyr::select(
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
    Interpretation
  )
print(publication_table)
write.csv(publication_table, "CE_AB_analysis/ab_triage_publication_table.csv", row.names = FALSE)

absolute_table <- results %>%
  mutate(
    `Monitoring interval` = freq_label,
    `Annual cost: PCR-only (US$)` = round(cost_pcronly, 2),
    `Suppression-years: PCR-only` = round(eff_pcronly, 4),
    `Mean time unsuppressed: PCR-only (years)` = round(unsuppressed_pcronly, 4),
    `Mortality risk: PCR-only` = round(mort_pcronly, 4),
    `DALYs: PCR-only` = round(DALY_pcronly, 4),
    `Annual cost: Ab-triage (US$)` = round(cost_triage, 2),
    `Suppression-years: Ab-triage` = round(eff_triage, 4),
    `Mean time unsuppressed: Ab-triage (years)` = round(unsuppressed_triage, 4),
    `Mortality risk: Ab-triage` = round(mort_triage, 4),
    `DALYs: Ab-triage` = round(DALY_triage, 4)
  ) %>%
  dplyr::select(
    `Monitoring interval`, `Annual cost: PCR-only (US$)`, `Suppression-years: PCR-only`,
    `Mean time unsuppressed: PCR-only (years)`, `Mortality risk: PCR-only`, `DALYs: PCR-only`,
    `Annual cost: Ab-triage (US$)`, `Suppression-years: Ab-triage`,
    `Mean time unsuppressed: Ab-triage (years)`, `Mortality risk: Ab-triage`, `DALYs: Ab-triage`
  )
print(absolute_table)
write.csv(absolute_table, "CE_AB_analysis/absolute_table_publication_table.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 6) Point-of-Care (POC) Comparison & Deterministic Panel Figure
# ------------------------------------------------------------------------------
results_poc <- scenario_grid %>%
  rowwise() %>%
  mutate(out = list(calc_per_person(
    frequency_yrs = frequency_yrs, 
    params_list = params, 
    annual_non_suppression = annual_non_suppression, 
    cost_visit_ab = 0
  ))) %>%
  select(-any_of("strategy")) %>% # Drop duplicate column before unnesting
  unnest(out) %>%
  mutate(
    freq_label = case_when(
      frequency_yrs == 1 ~ "Annual",
      frequency_yrs == 0.5 ~ "Biannual",
      frequency_yrs == 0.25 ~ "Quarterly",
      frequency_yrs == 0.125 ~ "6-weekly",
      abs(frequency_yrs - 1/12) < 1e-6 ~ "Monthly",
      TRUE ~ paste0(round(1/frequency_yrs, 1), "x/yr")
    ),
    freq_label = factor(freq_label, levels = c("6-weekly", "Quarterly", "Biannual", "Annual", "Monthly")),
    scenario = "Point-of-care"
  )

results_facility <- results %>% mutate(scenario = "Facility-based")

results_compare <- bind_rows(results_facility, results_poc) %>%
  mutate(NMB_500 = 500 * delta_DALY - delta_cost, prevalence_label = paste0(annual_non_suppression * 100, "% prevalence"))

p_poc_nmb <- ggplot(results_compare, aes(x = freq_label, y = NMB_500, group = scenario, color = scenario)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  facet_wrap(~ prevalence_label) +
  scale_color_manual(values = c("Facility-based" = "#d95f02", "Point-of-care" = "#1b9e77")) +
  labs(x = "Monitoring interval", y = "Net monetary benefit (US$)", color = "", title = "Economic impact of point-of-care antibody monitoring") +
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"), legend.position = "bottom")

print(p_poc_nmb)
ggsave("CE_AB_analysis/p_poc_nmb.png", p_poc_nmb, width = 10, height = 8)

# Deterministic 4-Panel Figure
deterministic_plot_data <- results %>%
  mutate(prevalence_label = paste0(round(annual_non_suppression * 100), "%"), NMB_500 = 500 * delta_DALY - delta_cost) %>%
  dplyr::select(freq_label, prevalence_label, delta_cost, delta_DALY, ICER_DALY, NMB_500) %>%
  pivot_longer(cols = c(delta_cost, delta_DALY, ICER_DALY, NMB_500), names_to = "Outcome", values_to = "Value") %>%
  mutate(Outcome = recode(Outcome, delta_cost = "Incremental cost (US$)", delta_DALY = "DALYs averted", ICER_DALY = "ICER (US$/DALY averted)", NMB_500 = "Net monetary benefit (US$)"))

plot_cost <- ggplot(filter(deterministic_plot_data, Outcome == "Incremental cost (US$)"), aes(freq_label, Value, group = prevalence_label, color = prevalence_label)) +
  geom_line(size = 1.2) + geom_point(size = 3) + labs(x = "Monitoring interval", y = "US$", color = "Prevalence") + theme_minimal(base_size = 14)

plot_daly <- ggplot(filter(deterministic_plot_data, Outcome == "DALYs averted"), aes(freq_label, Value, group = prevalence_label, color = prevalence_label)) +
  geom_line(size = 1.2) + geom_point(size = 3) + labs(x = "Monitoring interval", y = "DALYs", color = "Prevalence") + theme_minimal(base_size = 14)

plot_icer <- ggplot(filter(deterministic_plot_data, Outcome == "ICER (US$/DALY averted)"), aes(freq_label, Value, group = prevalence_label, color = prevalence_label)) +
  geom_line(size = 1.2) + geom_point(size = 3) + labs(x = "Monitoring interval", y = "US$/DALY", color = "Prevalence") + theme_minimal(base_size = 14)

plot_nmb <- ggplot(filter(deterministic_plot_data, Outcome == "Net monetary benefit (US$)"), aes(freq_label, Value, group = prevalence_label, color = prevalence_label)) +
  geom_line(size = 1.2) + geom_point(size = 3) + labs(x = "Monitoring interval", y = "US$", color = "Prevalence") + theme_minimal(base_size = 14)

final_deterministic_figure <- (plot_cost + plot_daly) / (plot_icer + plot_nmb) + plot_annotation(tag_levels = "A")
print(final_deterministic_figure)
ggsave("CE_AB_analysis/deterministic_headline_figure.png", final_deterministic_figure, width = 12, height = 8, dpi = 600)

# ------------------------------------------------------------------------------
# 7) Value-of-Information (VOI) Decision Contours
# ------------------------------------------------------------------------------
voi_grid <- expand.grid(
  ab_sens = seq(0.50, 0.95, length.out = 30),
  cost_ab = seq(0.50, 10.0, length.out = 30),
  frequency_yrs = 0.25
) %>% as_tibble()

voi_results <- voi_grid %>%
  rowwise() %>%
  mutate(out = list(calc_per_person(frequency_yrs = frequency_yrs, strategy = "AB_triage", ab_sens = ab_sens, cost_ab = cost_ab))) %>%
  unnest(out) %>%
  mutate(
    NMB_500 = 500 * delta_DALY - delta_cost,
    CE_Category = case_when(
      delta_cost < 0 & delta_DALY > 0 ~ "Dominant",
      NMB_500 > 0 ~ "Cost-effective",
      delta_DALY < 0 ~ "Health-harming",
      TRUE ~ "Not cost-effective"
    )
  )

fig_voi_contour <- ggplot(voi_results, aes(x = ab_sens, y = cost_ab, fill = CE_Category)) +
  geom_tile(alpha = 0.8) +
  scale_fill_manual(values = c("Dominant" = "#1b9e77", "Cost-effective" = "#7fc97f", "Not cost-effective" = "#fdc086", "Health-harming" = "#d95f02")) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = dollar_format()) +
  labs(x = "Antibody Test Sensitivity", y = "Unit Cost of Antibody Test (USD)", fill = "Economic Verdict", title = "Decision Boundaries for Quarterly Antibody Triage (WTP = $500/DALY)") +
  theme_minimal(base_size = 14)

print(fig_voi_contour)
ggsave("CE_AB_analysis/voi_decision_boundaries.png", fig_voi_contour, width = 9, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 8) One-Way Tornado Sensitivity Analysis
# ------------------------------------------------------------------------------
tornado_inputs <- tibble(
  param = c("ab_sensitivity", "ab_specificity", "cost_ab", "cost_pcr", "mean_gain_delay_reduction", "prop_resistance", "p_resupp_after_resistance_detect", "p_resupp_after_behav_detect", "annual_non_suppression", "mort_rate_suppressed", "mort_rate_unsuppressed"),
  base  = c(params$ab_sensitivity, params$ab_specificity, params$cost_ab, params$cost_pcr, params$mean_gain_delay_reduction, params$prop_resistance, params$p_resupp_after_resistance_detect, params$p_resupp_after_behav_detect, params$annual_non_suppression, params$mort_rate_suppressed, params$mort_rate_unsuppressed),
  low   = c(0.3, 0.85, 1.0, 20, 0.05, 0.1, 0.05, 0.2, 0.01, 0.005, 0.03),
  high  = c(0.9, 0.99, 7.0, 100, 0.5, 0.7, 0.5, 0.95, 0.30, 0.03, 0.10)
)

tornado_all <- freq_scenarios %>%
  rowwise() %>%
  mutate(
    tornado = list({
      base_freq <- frequency_yrs
      tornado_inputs %>%
        rowwise() %>%
        mutate(
          out_low = list({ pmod <- params; pmod[[param]] <- low; calc_per_person(base_freq, params_list = pmod) }),
          out_high = list({ pmod <- params; pmod[[param]] <- high; calc_per_person(base_freq, params_list = pmod) })
        ) %>%
        unnest(cols = c(out_low, out_high), names_sep = "_") %>%
        mutate(
          NMB_low  = 500 * out_low_delta_DALY  - out_low_delta_cost,
          NMB_high = 500 * out_high_delta_DALY - out_high_delta_cost,
          NMB_base = { base <- calc_per_person(base_freq, params_list = params); 500 * base$delta_DALY - base$delta_cost },
          min_change = pmin(NMB_low - NMB_base, NMB_high - NMB_base),
          max_change = pmax(NMB_low - NMB_base, NMB_high - NMB_base)
        )
    })
  ) %>%
  unnest(tornado)

tornado_plot_all <- ggplot(tornado_all) +
  geom_linerange(aes(x = param, ymin = min_change, ymax = max_change), size = 5) +
  coord_flip() +
  facet_wrap(~ forcats::fct_relevel(freq_label, "6-weekly", "Quarterly", "Biannual", "Annual"), scales = "free_x") +
  labs(x = "", y = "Change in NMB (USD, WTP = 500 per DALY averted)", title = "One-way sensitivity analysis across monitoring frequencies") +
  theme_minimal()

print(tornado_plot_all)
ggsave("CE_AB_analysis/tornado_all_frequencies.png", tornado_plot_all, width = 12, height = 8)

# ------------------------------------------------------------------------------
# 9) Probabilistic Sensitivity Analysis (PSA)
# ------------------------------------------------------------------------------
n_sims <- 2000
draws <- tibble(
  ab_sensitivity = rbeta(n_sims, 60, 40),
  ab_specificity = rbeta(n_sims, 95, 5),
  cost_ab = rgamma(n_sims, shape = 4, scale = 0.5),
  cost_pcr = rgamma(n_sims, shape = 40, scale = 1.0),
  prop_res = rbeta(n_sims, 40, 60),
  p_resup_res = rbeta(n_sims, 20, 80),
  p_resup_beh = rbeta(n_sims, 60, 40),
  annual_non_suppression = rbeta(n_sims, 10, 90),
  mort_rate_suppressed = rbeta(n_sims, 2, 198),
  mort_rate_unsuppressed = rbeta(n_sims, 6, 94)
)

WTP_threshold <- 500
psa_all <- freq_scenarios %>%
  rowwise() %>%
  mutate(
    psa = list({
      psa_freq <- frequency_yrs
      draws %>%
        mutate(sim = row_number()) %>%
        rowwise() %>%
        mutate(out = list(calc_per_person(
          frequency_yrs = psa_freq,
          strategy = "AB_triage",
          annual_non_suppression = annual_non_suppression,
          ab_sens = ab_sensitivity,
          ab_spec = ab_specificity,
          cost_ab = cost_ab,
          cost_pcr = cost_pcr,
          prop_res = prop_res,
          p_resup_res = p_resup_res,
          p_resup_beh = p_resup_beh,
          mort_rate_suppressed = mort_rate_suppressed,
          mort_rate_unsuppressed = mort_rate_unsuppressed
        ))) %>%
        unnest(out) %>%
        transmute(delta_cost, delta_eff, delta_LY, delta_DALY)
    })
  ) %>%
  unnest(psa) %>%
  mutate(
    delta_days = delta_DALY * 365,
    CE_status = ifelse(delta_cost <= WTP_threshold * delta_DALY, "Cost-effective", "Not cost-effective"),
    panel_label = case_when(
      freq_label == "6-weekly" ~ "A. Every 6 weeks",
      freq_label == "Quarterly" ~ "B. Every 3 months",
      freq_label == "Biannual" ~ "C. Every 6 months",
      freq_label == "Annual" ~ "D. Every 12 months"
    )
  )

ce_plane_all <- ggplot(psa_all, aes(x = delta_days, y = delta_cost, color = CE_status)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_abline(intercept = 0, slope = WTP_threshold / 365, linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.4) +
  facet_wrap(~ panel_label) +
  scale_color_manual(values = c("Cost-effective" = "#1b9e77", "Not cost-effective" = "#d95f02")) +
  labs(x = "Health gain (DALYs averted, days)", y = "Additional cost (US$)", color = "Simulation result") +
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 15, face = "bold"), legend.position = "bottom")

print(ce_plane_all)
ggsave("CE_AB_analysis/ce_plane_all_frequencies.png", ce_plane_all, width = 12, height = 8)

# CEAC
ceac_all <- psa_all %>%
  crossing(WTP = seq(0, 1000, by = 20)) %>%
  mutate(NMB = WTP * delta_DALY - delta_cost) %>%
  group_by(freq_label, WTP) %>%
  summarize(pr_CE = mean(NMB > 0), .groups = "drop") %>%
  mutate(
    panel_label = case_when(
      freq_label == "6-weekly" ~ "A. Every 6 weeks",
      freq_label == "Quarterly" ~ "B. Every 3 months",
      freq_label == "Biannual" ~ "C. Every 6 months",
      freq_label == "Annual" ~ "D. Every 12 months"
    ),
    pr_CE_percent = pr_CE * 100
  )

ceac_plot_all <- ggplot(ceac_all, aes(x = WTP, y = pr_CE_percent)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ factor(panel_label, levels = c("A. Every 6 weeks", "B. Every 3 months", "C. Every 6 months", "D. Every 12 months"))) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(x = "Willingness-to-pay threshold (US$ per DALY averted)", y = "Probability the strategy is cost-effective") +
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 15, face = "bold"))

print(ceac_plot_all)
ggsave("CE_AB_analysis/ceac_all_frequencies.png", ceac_plot_all, width = 12, height = 8)

message("Analysis completed successfully. All figures and tables exported to 'CE_AB_analysis/'.")
