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

baseline_ODn_table <- full_dataset %>%
  group_by(subject_label_blinded) %>%
  mutate(n_visit = 1:length(subject_label_blinded)) %>%
  ungroup() %>%
  filter(n_visit == 1)
mean(c(baseline_ODn_table$sedia_ODn, baseline_ODn_data$ODn)) ## baseline mean
sd(c(baseline_ODn_table$sedia_ODn, baseline_ODn_data$ODn))

# model_eddi <- nlme::lme(sedia_ODn ~ poly(days_since_eddi, 2), 
#                   data = sedia_eddi_data, 
#                   random = ~ days_since_eddi|subject_label_blinded,
#                   na.action = na.omit)
# summary(model_eddi)
# ploy_model <- nlme::lme(sedia_ODn ~ poly(years_since_tx_start, 2), 
#                         data = full_dataset %>%
#                           filter(Group == 'early suppressions'), 
#                         random = ~years_since_tx_start|subject_label_blinded,
#                         na.action = na.omit)
simple_model <- lm(sedia_ODn ~ years_since_tx_start, 
                   data = full_dataset %>%
                     mutate(sedia_ODn = log(sedia_ODn)))

# Get the coefficients
start_values <- coef(simple_model)

# Use the coefficients as starting values
poly_model <- glmmTMB(
  sedia_ODn ~ years_since_tx_start + (1 | subject_label_blinded),
  family = gaussian(link = "log"),
  data = full_dataset,
  start = list(beta = c(0, 1))
)
coef_estimates <- coefs[, "Estimate"]
coef_se <- coefs[, "Std. Error"]
summary(poly_model)
# summary(poly_model)$tTable[,2][[1]]
# sd_fixed <- c(summary(poly_model)$tTable[,2][[1]], 
#               summary(poly_model)$tTable[,2][[2]],
#               summary(poly_model)$tTable[,2][[3]]) * sqrt(length(unique(sedia_eddi_data$subject_label_blinded)))

#################################################################################
####individual decay rates
###################################################################################

# Set seed for reproducibility
set.seed(123)

# Number of individuals
n_individuals <- 10000

# Time points (6-month intervals over 10 years)
time_points <- seq(0, 10, by = 0.5)

# Define the distributions for baseline
baseline_mean <- 3.47
baseline_sd <- 1.55
baselines <- truncnorm::rtruncnorm(n_individuals, mean = baseline_mean, sd = baseline_sd, a = 0.03, b = 7.4)

# Define a two-degree polynomial to generate decay parameters
generate_decay_rate <- function(t) {
  # browser()
  # Two-degree polynomial: a*t^2 + b*t + c
  a <- rnorm(1, mean = coef_estimates[[1]], sd = coef_se[[1]])#1.818730
  b <- rnorm(1, mean = coef_estimates[[2]], sd = coef_se[[2]])#-8.97018
  decay_rate <- a + b * t
  return(c(pmax(decay_rate, 0.05), a, b))  # Ensure decay rates are non-negative
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
# decay_data <- data.frame(time = rep(time_points, n_individuals),
#                          individual = rep(1:n_individuals, each = length(time_points)),
#                          value = unlist(lapply(1:n_individuals, function(i) {
#                            exp_decay(time_points, baselines[i], decay_rates[1:length(n_individuals), i])
#                          })))

# Define measurement noise as a function of ODn
sigma_0 <- -0.01469  # Example intercept for noise model
slope_sigma <- 0.14513  # Example slope

# Function to compute noise SD based on ODn value
compute_noise_sd <- function(odn) {
  pmax(sigma_0 + slope_sigma * odn, 0.01)  # Ensure non-negative SD
}

# Generate decay data with heteroskedastic noise
decay_data <- data.frame(time = rep(time_points, n_individuals),
                         individual = rep(1:n_individuals, each = length(time_points)),
                         value = unlist(lapply(1:n_individuals, function(i) {
                           expected_odn <- exp_decay(time_points, baselines[i], decay_rates[1:length(time_points), i])
                           noise_sd <- sapply(expected_odn, compute_noise_sd)
                           noisy_odn <- pmax(rnorm(length(expected_odn), mean = expected_odn, sd = noise_sd), 0.001)
                           return(noisy_odn)
                         })))

model_parameters <- bind_cols(id = 1:n_individuals,
                              baseline = as.data.frame(baselines), 
                              remove_rownames(data.frame(t(data.frame(decay_rates[c(22,23,24),])))) %>%
                                dplyr::select(a = X1, b = X2, c = X3))

saveRDS(model_parameters, 'data/Exponential.rds')

x <- full_dataset %>%
  group_by(subject_label_blinded) %>%
  arrange(subject_label_blinded, years_since_tx_start) %>%
  mutate(visits = 1:length(subject_label_blinded),
         max_visits = max(visits)) %>%
  filter(Group == 'early suppressions') %>%
  distinct(subject_label_blinded, .keep_all = T)
x2 = cephia_samples %>%
  distinct(subject_label_blinded, .keep_all = T)

dataset_ggplot <- decay_data %>%
  filter(individual %in% sample(length(unique(decay_data$individual)), 100, replace = F)) 
set.seed(11)
ggplot_plots <- ggpubr::ggarrange( 
  ggplot(full_dataset %>%
           filter(Group == 'early suppressions'), 
         aes(x = years_since_tx_start, y = sedia_ODn, group = subject_label_blinded, color = subject_label_blinded)) +
    geom_line(size = 1.5) +
    # geom_line(alpha = 0.5) +
    labs(#title = "Exponential Decay Curves for Individuals",
      x = "Time since ART start (years)",
      y = "ODn Value") +
    # theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      legend.position = "none"
    ),
  
  
  ggplot(decay_data %>%
           mutate(flag = as.logical(ifelse(individual %in% sample(n_individuals, 50, replace = F), 1, 0))), 
         aes(x = time, y = value, group = individual, color = flag)) +
    geom_line(alpha = 0.5, linewidth = 1.5) +
    gghighlight(flag, use_direct_label = FALSE, unhighlighted_colour = "grey70") +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "black")) +
    labs(#title = "Exponential Decay Curves for Individuals",
      x = "Time since ART start (years)",
      y = "Estimated ODn Value")  +
    scale_y_continuous(limits = c(0, max(dataset_ggplot$value))) +
    scale_x_continuous(limits = c(0, max(dataset_ggplot$time))) +
    # theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      legend.position = "none"
    ),
  labels = c("A", "B"),
  ncol = 1, nrow = 2
)

# ggplot(decay_data %>%
#          mutate(flag = as.logical(ifelse(individual %in% sample(n_individuals, 50, replace = F), 1, 0))), 
#        aes(x = time, y = value, group = individual, color = flag)) +
#   geom_line(alpha = 0.5, linewidth = 1.5) +
#   gghighlight(flag, use_direct_label = FALSE, unhighlighted_colour = "grey70") +
#   scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "black")) +
#   labs(title = "Exponential Decay Curves for Individuals",
#        x = "Time since ART start (years)",
#        y = "ODn Value") +
#   theme_minimal() +
#   theme(legend.position="none")

jpeg('other_figures/simulated_plot - exponential.jpeg', units = "in", width = 9, height = 9, res = 300)
ggplot_plots
dev.off()

###########################################################
##getting data to use for testing - CEPHIA EDCTP data
###########################################################
dt02 <- read_csv("data/20180410-EP-LAgSedia-Generic.csv") %>%
  filter(on_treatment == T) %>%
  dplyr::select(subject_label_blinded, visit_date, test_date, art_initiation_date, 
                days_since_eddi, viral_load, sedia_ODn = result...72) %>% 
  arrange(subject_label_blinded, test_date) %>%
  # filter(as.character(art_initiation_date) != '') %>%
  # filter(visit_date >= art_initiation_date) %>%
  mutate(years_since_ART_start = as.numeric(visit_date-art_initiation_date)/365.25) %>%
  arrange(subject_label_blinded, visit_date) %>%
  group_by(subject_label_blinded) %>%
  mutate(x = 1:length(subject_label_blinded)) %>%
  mutate(flag = max(x)) %>%
  ungroup() %>%
  filter(flag >2) %>%
  group_by(subject_label_blinded) %>%
  mutate(x = 1:length(subject_label_blinded),
         flag = max(x)) %>%
  mutate(peak = ifelse(x==flag, 1, NA)) %>%
  ungroup() %>%
  mutate(strata = 'early suppressions - cephia') %>%
  dplyr::select(subject_label_blinded, strata, years_since_ART_start,sedia_ODn, peak)
# write.csv(x, 'output_table/visits_during_ART.csv')
# dt03 <- read_csv('output_table/visits_during_ART - edited.csv') %>%
#   filter(selected_visits == 1) %>%
#   mutate(years_since_eddi = days_since_eddi/365.25) %>%
#   dplyr::select(subject_label_blinded, years_since_eddi, sedia_ODn, selected_visits, peak)
dt04 <- full_dataset %>%
  filter(Group == 'suppression failures') %>%
  group_by(subject_label_blinded) %>%
  mutate(x = 1:length(subject_label_blinded),
         flag = max(x)) %>%
  mutate(peak = ifelse(x==flag, 1, NA)) %>%
  ungroup() %>%
  filter(flag >2) %>%
  mutate(years_since_ART_start = `days since tx start`/365.25) %>%
  dplyr::select(subject_label_blinded, strata = Group, years_since_ART_start, sedia_ODn, peak)
test_data <- bind_rows(
  dt02 %>%
    mutate(subject_label_blinded = as.character(subject_label_blinded)), dt04
) %>%
  mutate(id = as.numeric(as.factor(subject_label_blinded))) %>%
  dplyr::select(record_id = id, time_vec = years_since_ART_start, 
                ODn_vec = sedia_ODn, peak_visit = peak)

# best_model_choice <- function(test_data, param_sim_data) {
#   # browser()
#   time_vec <- test_data[,2]
#   ODn_vec <- test_data[,3]
#   dt <- model_parameters
#   a <- model_parameters$a
#   b <- model_parameters$b
#   c <- model_parameters$c
#   for (i in 1:length(ODn_vec[[1]])) {
#     t <- time_vec[[i,1]]
#     y <- ODn_vec[[i,1]]
#     dt[[paste0("square_error", i)]] <- (y - pmax(baselines * exp(-(a * t^2 + b * t + c) * t), .001))^2
#   }
#   # browser()
#   dt1 <- dt[!is.infinite(rowSums(dt)),]
#   # stopifnot(length(dt1$baselines) == 0)
#   dt2 <- dt1 %>%
#     rowwise() %>%
#     mutate(
#       rmse = sqrt(mean(c_across(starts_with("square_error"))))
#     ) %>%
#     ungroup() %>%
#     mutate(
#       baselines = as.numeric(baselines),  # <-- Add this line
#       min_slope = min(rmse)
#     ) %>%
#     filter(rmse == min(rmse)) %>%
#     dplyr::select(id, baselines, a, b, c, rmse)
#   return(dt2)
# }

# results <- c()
# ids <- unique(test_data$record_id)
# for (i in 1:length(unique(test_data$record_id))) {
#   # print(ids[i])
#   best_model_parameters <- best_model_choice(test_data = test_data %>% 
#                                                filter(is.na(peak_visit)) %>%
#                                                filter(record_id == ids[i]), 
#                                              param_sim_data = model_parameters)
#   results <- rbind(results, best_model_parameters)
# }
best_model_choice <- function(test_data, param_sim_data, baselines) {
  time_vec <- test_data[, 2]
  ODn_vec <- test_data[, 3]
  
  dt <- param_sim_data
  a <- dt$a
  b <- dt$b
  c <- dt$c
  
  # Loop over each timepoint
  for (i in 1:nrow(test_data)) {
    t <- time_vec[[i, 1]]
    y <- ODn_vec[[i, 1]]
    pred <- pmax(baselines * exp(-(a * t^2 + b * t + c) * t), 0.001)
    dt[[paste0("squared_error", i)]] <- (y - pred)^2
    dt[[paste0("abs_error", i)]] <- abs(y - pred)
  }
  
  # Filter out rows with non-finite values
  dt_clean <- dt[!is.infinite(rowSums(dt %>% dplyr::select(starts_with("squared_error")))), ]
  
  # Compute RMSE and MAE
  dt_scores <- dt_clean %>%
    rowwise() %>%
    mutate(
      rmse = sqrt(mean(c_across(starts_with("squared_error")))),
      mae = mean(c_across(starts_with("abs_error"))),
      complexity = abs(a) + abs(b) + abs(c)
    ) %>%
    ungroup()
  
  # Step 1: Filter for lowest RMSE
  min_rmse <- min(dt_scores$rmse)
  dt_min_rmse <- dt_scores %>% filter(rmse == min_rmse)
  
  # Step 2: If ties, use MAE
  if (nrow(dt_min_rmse) > 1) {
    min_mae <- min(dt_min_rmse$mae)
    dt_min_rmse <- dt_min_rmse %>% filter(mae == min_mae)
  }
  
  # Step 3: If still ties, use lowest model complexity
  if (nrow(dt_min_rmse) > 1) {
    min_complexity <- min(dt_min_rmse$complexity)
    dt_min_rmse <- dt_min_rmse %>% filter(complexity == min_complexity)
  }
  
  # Step 4: Still ties? Pick one randomly
  if (nrow(dt_min_rmse) > 1) {
    set.seed(11)
    dt_min_rmse <- dt_min_rmse[sample(1:nrow(dt_min_rmse), 1), ]
  }
  
  # Final selection
  return(dt_min_rmse %>% dplyr::select(id, baselines, a, b, c, rmse, mae))
}


results <- list()
ids <- unique(test_data$record_id)

for (i in seq_along(ids)) {
  individual_data <- test_data %>%
    filter(is.na(peak_visit)) %>%
    filter(record_id == ids[i])
  
  best_model_parameters <- best_model_choice(
    test_data = individual_data,
    param_sim_data = model_parameters,
    baselines = model_parameters$baselines
  )
  
  results[[i]] <- best_model_parameters
}

results_df <- do.call(rbind, results)

results <- results_df %>%
  mutate(record_id = ids)
head(results)
saveRDS(results, 'results_100k.rds')
results <- readRDS('results_100k.rds')
###########################################################
###sigma ODn
############################################################
controls_blinded_sedia <-  read_csv("data/20180410-EP-LAgSedia-Generic.csv") %>%
  mutate(sedia_ODn = `result...15`) %>%
  filter(visit_id %in% c(21773, 21783, 21785)) %>%
  mutate(specimen = case_when(
    visit_id == 21773 ~ "BC-1",
    visit_id == 21783 ~ "BC-2",
    visit_id == 21785 ~ "BC-3"
  )) %>%
  select(specimen, test_date,
         specimen_blinded_label = specimen_label,
         sedia_final_ODn = `result...15`,
         sedia_final_ODn_method = method
  )
sedia_distribution_blinded_controls <- controls_blinded_sedia %>%
  group_by(specimen) %>%
  summarise(
    `mean Sedia ODn` = round(mean(sedia_final_ODn), 3),
    `sigma Sedia ODn` = round(sd(sedia_final_ODn), 3)
  ) %>%
  mutate(`sigma over mean ODn` = `sigma Sedia ODn` / `mean Sedia ODn`)

sd_sedia_ODn_vs_ODn_bc <- sedia_distribution_blinded_controls %>%
  ggplot(aes(x = `mean Sedia ODn`, y = `sigma Sedia ODn`)) +
  geom_point() +
  expand_limits(x = 0, y = 0) +
  geom_smooth(method = lm, size = 1.5, se = FALSE)

sd_sedia_ODn_vs_ODn_bc

noise <- summary(glm(sd_Sedia_ODn ~ mean_Sedia_ODn, data = sedia_distribution_blinded_controls %>%
                       mutate(sd_Sedia_ODn = `sigma Sedia ODn`, mean_Sedia_ODn = `mean Sedia ODn`)))
######################################################################
#How about a new time point?
######################################################################
compare_value_with_others <- function(data_set, t, y_ODn, sigma_ODn, sigma_y_ODn, id) {
  # browser()
  t <- t #/ 365.25
  y <- y_ODn
  y_hat <- (data_set %>%
              mutate(y_hat = pmax(baselines * exp(-(a * t^2 + b * t + c) * t), .001)))$y_hat
  sigma_pooled <- (sigma_y_ODn^2 + sigma_ODn^2)^.5
  z_test <- (y - y_hat) / sigma_pooled
  set.seed(11)
  Z <- Normal(0, 1) # make a standard normal r.v.
  p_value <- 1 - pnorm(z_test) #2 * (1 - cdf(Z, abs(z_test)))
  results <- cbind(id = data_set$id, value = y, z_stat = z_test, p_value = p_value, record_id = id) # value = y,
  return(results)
}
test_data_last_visit <- test_data %>%
  filter(!is.na(peak_visit))
results1 <- c()
for (i in 1:length(results$record_id)) {
  results1 <- rbind(results1, compare_value_with_others(data_set = results %>% filter(record_id == results$record_id[i]),#best_model_parameters, 
                                                        t = test_data_last_visit$time_vec[test_data_last_visit$record_id == results$record_id[i]], #Must be in days
                                                        y_ODn = test_data_last_visit$ODn_vec[test_data_last_visit$record_id == results$record_id[i]], 
                                                        sigma_ODn = sd(test_data$ODn_vec[test_data$record_id == results$record_id[i] & is.na(test_data$peak_visit)]),
                                                        sigma_y_ODn = noise$coefficients[1,1] + noise$coefficients[2,1] * (test_data_last_visit$time_vec[test_data_last_visit$record_id == results$record_id[i]]), #sd((full_dataset %>% filter(Group == 'early suppressions'))$sedia_ODn), #sd(dt$ODn),
                                                        id = results$record_id[i]
  ))
}
dt05 <- as.data.frame(results1) %>%
  left_join(bind_rows(
    dt02 %>%
      mutate(subject_label_blinded = as.character(subject_label_blinded)), dt04
  ) %>%
    mutate(record_id = as.numeric(as.factor(subject_label_blinded))), by = 'record_id') %>%
  dplyr::select(names(as.data.frame(results1)), strata) %>%
  distinct(record_id, .keep_all = T) %>%
  mutate(y_hat_status = as.factor(ifelse(p_value < 0.05, 1, 2)),
         y = as.factor(ifelse(strata == 'early suppressions', 1, 2)) )

dt06 <- dt05 %>%
  dplyr::select(strata, y_hat_status) %>%
  tbl_summary(by = strata
  )
dt06
dt07 <- dt05 %>%
  mutate(`significance 99%` = ifelse(as.numeric(z_stat) > qnorm(0.99), TRUE, FALSE),
         `significance 98.5%` = ifelse(as.numeric(z_stat) > qnorm(0.985), TRUE, FALSE),
         `significance 98%` = ifelse(as.numeric(z_stat) > qnorm(0.98), TRUE, FALSE),
         `significance 97.5%` = ifelse(as.numeric(z_stat) > qnorm(0.975), TRUE, FALSE),
         `significance 95%` = ifelse(as.numeric(z_stat) > qnorm(0.95), TRUE, FALSE),
         `P-value` = as.numeric(p_value)) %>%
  filter(strata == 'early suppressions - cephia')
dt08 <- dt05 %>%
  mutate(`significance 99%` = ifelse(as.numeric(z_stat) > qnorm(0.99), TRUE, FALSE),
         `significance 98.5%` = ifelse(as.numeric(z_stat) > qnorm(0.985), TRUE, FALSE),
         `significance 98%` = ifelse(as.numeric(z_stat) > qnorm(0.98), TRUE, FALSE),
         `significance 97.5%` = ifelse(as.numeric(z_stat) > qnorm(0.975), TRUE, FALSE),
         `significance 95%` = ifelse(as.numeric(z_stat) > qnorm(0.95), TRUE, FALSE),
         `P-value` = as.numeric(p_value)) %>%
  filter(strata != 'early suppressions - cephia')

# Specificity
round(table(dt07$`significance 95%`)[[1]] / (table(dt07$`significance 95%`)[[1]] + table(dt07$`significance 95%`)[[2]]), 3)
specificity_value <- cbind(`Accuracy values` = c(
  round(table(dt07$`significance 95%`)[[1]] / (table(dt07$`significance 95%`)[[1]] + table(dt07$`significance 95%`)[[2]]), 3),
  round(table(dt07$`significance 97.5%`)[[1]] / (table(dt07$`significance 97.5%`)[[1]] + table(dt07$`significance 97.5%`)[[2]]), 3),
  round(table(dt07$`significance 98%`)[[1]] / (table(dt07$`significance 98%`)[[1]] + table(dt07$`significance 98%`)[[2]]), 3),
  round(table(dt07$`significance 98.5%`)[[1]] / (table(dt07$`significance 98.5%`)[[1]] + table(dt07$`significance 98.5%`)[[2]]), 3),
  round(table(dt07$`significance 99%`)[[1]] / (table(dt07$`significance 99%`)[[1]] + table(dt07$`significance 99%`)[[2]]), 3)
), level = c('Significance 95%', 'Significance 97.5%', 'Significance 98%', 'Significance 98.5%', 'Significance 99%'), Accuracy = rep(c('Specificity'), 5)
)
# Sensitivity
round(0 / (table(dt08$`significance 95%`)[[1]] + 0), 3) #table(dt08$`significance 95%`)[[1]] table(dt08$`significance 95%`)[[2]]
sensitivity_value <- cbind(`Accuracy values` = c(
  round(table(dt08$`significance 95%`)[[2]] / (table(dt08$`significance 95%`)[[1]] + table(dt08$`significance 95%`)[[2]]), 3),
  round(table(dt08$`significance 97.5%`)[[2]] / (table(dt08$`significance 97.5%`)[[1]] + table(dt08$`significance 97.5%`)[[2]]), 3),
  round(table(dt08$`significance 98%`)[[2]] / (table(dt08$`significance 98%`)[[1]] + table(dt08$`significance 98%`)[[2]]), 3),
  round(table(dt08$`significance 98.5%`)[[2]] / (table(dt08$`significance 98.5%`)[[1]] + table(dt08$`significance 98.5%`)[[2]]), 3),
  round(table(dt08$`significance 99%`)[[2]] / (table(dt08$`significance 99%`)[[1]] + table(dt08$`significance 99%`)[[2]]), 3)
), level = c('Significance 95%', 'Significance 97.5%', 'Significance 98%', 'Significance 98.5%', 'Significance 99%'), Accuracy = rep('Sensitivity', 5)
)
dt09 <- bind_rows(dt07, dt08)
accuracy_dataset <- c()
threshold <- seq(0,3, 0.2)
for (i in 1:length(threshold)) {
  dat_set1 <- dt09 %>%
    filter(strata != 'early suppressions') %>%
    mutate(x = ifelse(as.numeric(z_stat) > threshold[i], TRUE, FALSE))
  # if(length(table(dat_set1$x)) == 1 & names(table(dat_set1$x))[1] == 'FALSE') {
  #   sensitivity_value <- 1} else if(length(table(dat_set1$x)) == 1 & names(table(dat_set1$x))[1] == 'TRUE'){ #& names(dat_set1$x)[1] == 'TRUE'
  #     sensitivity_value <- 1
  #   }else{
  sensitivity_value <- round(table(dat_set1$x)[[2]] / (table(dat_set1$x)[[1]] + table(dat_set1$x)[[2]]), 3)
  # }
  
  
  dat_set2 <- dt09 %>%
    filter(strata == 'early suppressions - cephia') %>%
    mutate(x = ifelse(as.numeric(z_stat) > threshold[i], TRUE, FALSE))
  # if(length(table(dat_set2$x)) == 1 & names(table(dat_set2$x))[1] == 'FALSE') {
  #   specificity_value <- 1} else if(length(table(dat_set2$x)) == 1 & names(table(dat_set2$x))[1] == 'TRUE'){
  #     specificity_value <- 0
  #   }else{
  specificity_value <- round(table(dat_set2$x)[[1]] / (table(dat_set2$x)[[1]] + table(dat_set2$x)[[2]]), 3)
  # }
  dt <- c(`Z score` = threshold[i], value1 = sensitivity_value, value2 = specificity_value)
  accuracy_dataset <- rbind(accuracy_dataset, dt)
}

accuracy_graph <- data.frame(accuracy_dataset) %>%
  pivot_longer(cols = c('value1', 'value2'),
               names_to = 'accuracy',
               values_to = 'Accuracy.values') %>%
  mutate(Accuracy = ifelse(accuracy=='value1', 'Sensitivity', 'Specificity')) %>%
  ggplot(aes(x = Z.score, y = Accuracy.values, group = Accuracy)) +
  geom_line(aes(color = Accuracy), size = 1.5) + #stat = "identity", position = position_dodge()
  geom_point(aes(color = Accuracy), size = 3) +
  scale_color_manual(values=c('#999999','#E69F00'))  +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = c(seq(0, 3, 0.5)), labels = c(seq(0, 3, 0.5))) +
  ylab('') + xlab('Z value threshold') +
  labs(colour = NULL) +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "null")#,
    # axis.text.x = element_text(angle = 60, hjust = 1)
  )
accuracy_graph
