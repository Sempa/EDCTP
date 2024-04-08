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

# source('estimate_measurement_noise_dist.R')
#############################################
#' Estimating Sedia LAg Measurement noise
#############################################
cephia_pts <- read_csv("Sempa_final_pull_with_results.csv") %>%
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
  mutate(visits = 1:length(subject_label_blinded),
         unsupressed_visits = ifelse(viral_load >999, 1,0)) %>%
  mutate(visits = ifelse(unsupressed_visits==1, visits, NA),
         baseline_visit = ifelse(unsupressed_visits==1, max(visits, na.rm = T), 0)) %>%
  mutate(baseline_visit = ifelse((subject_label_blinded == 35329295 | subject_label_blinded == 52033420), 1, 
                                 ifelse(subject_label_blinded == 54382839, 5, 
                                        ifelse(subject_label_blinded == 64577680,13, baseline_visit)))) %>%
  mutate(baseline_visit = max(baseline_visit),
         visits = 1:length(subject_label_blinded)) %>%
  filter(visits >=baseline_visit) %>%
  mutate(`days since tx start` = `days since tx start` * -1,
         subject_label_blinded = as.character(subject_label_blinded),
         to_peak = NA) %>%
  select(subject_label_blinded, days_since_eddi, `days since tx start`, Sex = BiologicalSex, 
         Age = `Age at Draw`, test_date, sedia_ODn, viral_load, visits, Group, to_peak) # , baseline_visit 

full_dataset <- bind_rows(cephia_pts %>%
  mutate(cohort = 'cephia'),
  read_delim("data/full_africos_data_with_ODn.csv", 
             delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% #africos_pts <- 
    filter(exclude %in% c(1,2)) %>%
  mutate(days_since_eddi = NA,
         `days since tx start` = NA,
         Group = ifelse(exclude == 1, 'early suppressions', 
                        ifelse(exclude == 2, 'suppression failures', NA)),
         Group = ifelse(is.na(Group), 'suppression failures', Group),
         Age = ifelse(Age == '.', '31.6', Age)
         ) %>%
  select(subject_label_blinded = subjid, days_since_eddi, `days since tx start`, Sex, Age, 
         test_date = Date_Collected, sedia_ODn = ODn, viral_load = vl, 
         visits = Visit, Group, exclude, to_peak) %>%
  left_join(
readxl::read_excel("data/AFRICOS_Sempa Data Pull_24Jul23.xlsx", 
                           sheet = "Sheet1") %>%
  mutate(id = paste(Date_Collected, subjid, sep = '_')) %>%
  group_by(subjid) %>%
  mutate(switch_date = as.Date(as.numeric(ifelse(`ARV Start Date` == '.', '', `ARV Start Date`)), 
                               origin = "1900-01-01"),
         art_start = as.Date(min(as.numeric(ifelse(`ARV Start Date` == '.', '', `ARV Start Date`)), na.rm = T), 
                             origin = "1900-01-01")) %>%
  ungroup() %>%
  dplyr::select(subject_label_blinded = subjid, art_start) %>%
  distinct(subject_label_blinded, .keep_all = T), by = 'subject_label_blinded') %>%
  mutate(`days since tx start` = as.numeric(as.Date(test_date, origin = "1900-01-01") - art_start),
         cohort = 'africos',
         # subject_label_blinded = as.double(cur_group_id()),
         Age = as.numeric(Age),
         viral_load = as.numeric(viral_load))%>%
  select(subject_label_blinded, days_since_eddi, 
         `days since tx start`, Sex, Age, 
         test_date, sedia_ODn, viral_load, visits, Group, cohort, to_peak)
) %>% mutate(time_on_ART=`days since tx start`) %>%
  select(subject_label_blinded, days_since_eddi, 
         time_on_ART, Sex, Age, 
         test_date, sedia_ODn, viral_load, visits, Group, cohort, to_peak)

sedia_generic <- read_csv("data/20180410-EP-LAgSedia-Generic.csv") %>%
  mutate(sedia_ODn = `result...15`)

controls_blinded_sedia <- sedia_generic %>%
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

# sd_sedia_ODn_vs_ODn_bc

noise <- summary(glm(sd_Sedia_ODn ~ mean_Sedia_ODn, data = sedia_distribution_blinded_controls %>%
                       mutate(sd_Sedia_ODn = `sigma Sedia ODn`, mean_Sedia_ODn = `mean Sedia ODn`)))

######################################################
#' Analysis among VL suppressed at 100, 400, and 1000
######################################################

slopes_for_suppressed <- function(ODn_vl_data, threshold) {
  # browser()
  dir.create(paste("supp_", threshold, Sys.Date(), sep = "_"))
  folder_name <- paste("supp_", threshold, Sys.Date(), sep = "_")
  # browser()
  data_generated <- ODn_vl_data %>%
    select(
      subject_label_blinded, time_on_ART, test_date, sedia_ODn, viral_load # , art_initiation_date, aids_diagnosis_date,
      # art_interruption_date, art_resumption_date, treatment_naive,
      # on_treatment, first_treatment
    ) %>%
    arrange(subject_label_blinded, test_date) %>%
    # distinct(subject_label_blinded, test_date, .keep_all = T) %>%
    filter(!is.na(viral_load)) %>%
    group_by(subject_label_blinded) %>%
    mutate(flagvl = ifelse(viral_load <= threshold, 0, 1)) %>%
    mutate(suprressed_throughout_followup = ifelse(mean(flagvl) == 0, 1, 0)) %>%
    mutate(visits = 1:length(subject_label_blinded)) %>%
    mutate(n_visits = max(visits)) %>%
    filter(n_visits > 1) %>%
    ungroup()
  
  jpeg(paste(folder_name, "/", "Figure1", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)
  
  print(
    ggplot(
      data = data_generated %>%
        filter(suprressed_throughout_followup == 1) %>%
        group_by(subject_label_blinded),
      aes(x = time_on_ART, y = sedia_ODn)
    ) + # as.factor(subject_label_blinded)
      geom_line(aes(color = as.factor(subject_label_blinded)), size = 1) + # aes(color = subject_label_blinded, linetype = Group), group = Group
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
      ) +
      ylab("Sedia LAg ODn") +
      xlab("Time on ART (Days)") +
      ggtitle(paste("Sedia ODn during ART", threshold))
  )
  dev.off()
  # browser()
  model_data <- data.frame(
    id = NA, interceptlink_identity = NA, slope_link_identity = NA#,
    # intercept_link_log = NA, slope_link_log = NA
  )
  set.seed(11)
  dat <- data_generated %>%
    # filter(suprressed_throughout_followup == 1) %>%
    mutate(
      sedia_ODn_with_noise = sedia_ODn + rnorm(n = length(subject_label_blinded), mean = 0, sd = noise$coefficients[1, 1] + noise$coefficients[2, 1] * sedia_ODn),
      sedia_ODn_with_noise = ifelse(sedia_ODn_with_noise < 0, 0, sedia_ODn_with_noise),
      sedia_ODn_with_noise = ifelse(is.na(sedia_ODn_with_noise), sedia_ODn, sedia_ODn_with_noise)
    ) %>%
    filter(!is.na(time_on_ART)) %>%
    select(subject_label_blinded, test_date, time_on_ART, sedia_ODn, sedia_ODn_with_noise, viral_load)
  counter <- 0
  for (i in 1:length(unique(dat$subject_label_blinded))) {
    # if(i==13) browser()
    counter <- counter + 1
    model_data[counter, ] <- ODn_regression_function_supressed(
      dat = subset(dat, dat$subject_label_blinded == unique(dat$subject_label_blinded)[i])
    )
    # print(i)
  }
  model_data <- model_data %>%
    mutate(interceptlink_identity = as.numeric(interceptlink_identity),
           slope_link_identity = as.numeric(slope_link_identity))
  # browser()
  hist_plot1 <- hist(model_data$slope_link_identity, breaks = 30, freq = FALSE, ylab = '', xlab = 'Individual slopes', main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  area_hist1 <- sum(hist_plot1$counts* abs(hist_plot1$mids[1]-hist_plot1$mids[2]))
  n_values1 <- sum(hist_plot1$counts / area_hist1)
  set.seed(11)
  x <- rnorm(n_values1, mean = mean(model_data$slope_link_identity, na.rm = T), sd = sd(model_data$slope_link_identity))
  dx <- density(x)

  browser()#
  jpeg(paste(folder_name, "/", "Figure_combined", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 8, res = 300)
  slope_figure_1000 <- ggplot(
    data = model_data,
    aes(x = slope_link_identity)
  ) + # as.factor(subject_label_blinded)
    geom_histogram(color="black", fill="red", bins = 30, alpha = 0.9) +
    scale_x_continuous(limits = c(-0.006, 0.008), breaks = seq(-0.006, 0.007, by = 0.002), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
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
    ) +
    ylab("Count") +
    xlab("Slopes of individuals (ODn/Day)")
  
  slope_figure_1000
  
  dev.off()
  
  # browser()
  write_rds(model_data, 'data/model_data_art.rds')
  
  return(
    cbind(
      suppression = threshold,
      slope_mean_link_identity = mean(model_data$slope_link_identity),
      pvalue_link_identity_slope_diff_zero = t.test((model_data$slope_link_identity))$p.value#,
      # slope_mean_link_log = mean(model_data$slope_link_log),
      # pvalue_link_log_slope_diff_zero = t.test((model_data$slope_link_log))$p.value
    )
  )
}

ODn_regression_function_supressed <- function(dat, i=i) {
  # if(i==13)browser()
  model1 <- glm(
    formula = sedia_ODn_with_noise ~ time_on_ART,
    family = gaussian(link = "identity"),
    data = dat
  )
  
  # model2 <- glm(
  #   formula = sedia_ODn_with_noise ~ time_on_ART,
  #   family = gaussian(link = "log"),
  #   data = dat
  # )
  
  return(cbind(
    id = unique(dat$subject_label_blinded),
    intercept_link_identity = model1$coefficients[[1]],
    slope_link_identity = model1$coefficients[[2]]#,
    # intercept_link_log = model2$coefficients[[1]],
    # slope_link_log = model2$coefficients[[2]]
  ))
}


results <- c()
for (i in c(1000)) {# 100, 400, 
  x <- slopes_for_suppressed(ODn_vl_data = full_dataset %>% 
                               filter(Group == 'early suppressions') %>%
                               filter(time_on_ART >0 & viral_load < 1000), 
                             threshold = i)
  results <- rbind(results, x)
}
write.csv(results, "output_table/results_suppressed_art.csv") # results
set.seed(11);Z <- Normal(0, 1)
model_data <- read_rds('data/model_data_art.rds') %>%
  mutate(sigma_slope_identity = (sd(slope_link_identity)/(length(slope_link_identity) - 1))^0.5#,
         # sigma_slope_log = (sd(slope_link_log)/(length(slope_link_log) - 1))^0.5
         ) %>%
  group_by(id) %>%
  mutate(z_value_identity = (slope_link_identity - 0) / sigma_slope_identity,
         `P values identity` = 1 - cdf(Z, abs(z_value_identity)) + cdf(Z, -abs(z_value_identity)),
         # z_value_log = (slope_link_log - 0) / sigma_slope_log,
         # `P values log` = 1 - cdf(Z, abs(z_value_log)) + cdf(Z, -abs(z_value_log))
  )

################################################################
#' Analysis among early VL suppressed test straddling ART start
################################################################

slopes_for_suppressed <- function(ODn_vl_data, threshold) {
  # browser()
  dir.create(paste("supp_", threshold, Sys.Date(), sep = "_"))
  folder_name <- paste("supp_", threshold, Sys.Date(), sep = "_")
  # browser()
  data_generated <- ODn_vl_data %>%
    select(
      subject_label_blinded, time_on_ART, test_date, sedia_ODn, viral_load, ART_status # , art_initiation_date, aids_diagnosis_date,
      # art_interruption_date, art_resumption_date, treatment_naive,
      # on_treatment, first_treatment
    ) %>%
    arrange(subject_label_blinded, test_date) %>%
    # distinct(subject_label_blinded, test_date, .keep_all = T) %>%
    filter(!is.na(viral_load)) %>%
    group_by(subject_label_blinded) %>%
    mutate(flagvl = ifelse(viral_load <= threshold, 0, 1)) %>%
    mutate(suprressed_throughout_followup = ifelse(mean(flagvl) == 0, 1, 0)) %>%
    mutate(visits = 1:length(subject_label_blinded)) %>%
    mutate(n_visits = max(visits)) %>%
    filter(n_visits > 1) %>%
    ungroup()
  
  jpeg(paste(folder_name, "/", "Figure1_b4_ART", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)
  
  print(
    ggplot(
      data = data_generated %>%
        filter(suprressed_throughout_followup == 1) %>%
        group_by(subject_label_blinded),
      aes(x = time_on_ART, y = sedia_ODn)
    ) + # as.factor(subject_label_blinded)
      geom_line(aes(color = as.factor(subject_label_blinded)), size = 1) + # aes(color = subject_label_blinded, linetype = Group), group = Group
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
      ) +
      ylab("Sedia LAg ODn") +
      xlab("Time (Days)") +
      ggtitle(paste("Sedia ODn by ART status", threshold))
  )
  dev.off()
  # browser()
  model_data <- data.frame(
    id = NA, interceptlink_identity = NA, slope_link_identity = NA#,
    # intercept_link_log = NA, slope_link_log = NA
  )
  set.seed(11)
  dat <- data_generated %>%
    # filter(suprressed_throughout_followup == 1) %>%
    mutate(
      sedia_ODn_with_noise = sedia_ODn + rnorm(n = length(subject_label_blinded), mean = 0, sd = noise$coefficients[1, 1] + noise$coefficients[2, 1] * sedia_ODn),
      sedia_ODn_with_noise = ifelse(sedia_ODn_with_noise < 0, 0, sedia_ODn_with_noise),
      sedia_ODn_with_noise = ifelse(is.na(sedia_ODn_with_noise), sedia_ODn, sedia_ODn_with_noise)
    ) %>%
    filter(!is.na(time_on_ART)) %>%
    select(subject_label_blinded, test_date, time_on_ART, sedia_ODn, sedia_ODn_with_noise, viral_load, ART_status)
  counter <- 0
  # browser()
  for (i in 1:length(unique(dat$subject_label_blinded))) {
    # browser()
    counter <- counter + 1
    model_data[counter, ] <- ODn_regression_function_supressed(
      dat = subset(dat, dat$subject_label_blinded == unique(dat$subject_label_blinded)[i])
    )
    # print(i)
  }
  model_data <- model_data %>%
    mutate(interceptlink_identity = as.numeric(interceptlink_identity),
           slope_link_identity = as.numeric(slope_link_identity))
  # browser()
  hist_plot1 <- hist(model_data$slope_link_identity, breaks = 30, freq = FALSE, ylab = '', xlab = 'Individual slopes', main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  area_hist1 <- sum(hist_plot1$counts* abs(hist_plot1$mids[1]-hist_plot1$mids[2]))
  n_values1 <- sum(hist_plot1$counts / area_hist1)
  set.seed(11)
  x <- rnorm(n_values1, mean = mean(model_data$slope_link_identity, na.rm = T), sd = sd(model_data$slope_link_identity))
  dx <- density(x)
  
  browser()#
  jpeg(paste(folder_name, "/", "Figure_combined_b4_art", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 8, res = 300)
  slope_figure_1000 <- ggplot(
    data = model_data,
    aes(x = slope_link_identity)
  ) + # as.factor(subject_label_blinded)
    geom_histogram(color="black", fill="red", bins = 30, alpha = 0.9) +
    scale_x_continuous(limits = c(-0.006, 0.008), breaks = seq(-0.006, 0.007, by = 0.002), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
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
    ) +
    ylab("Count") +
    xlab("Slopes of individuals (ODn/Day)")
  
  slope_figure_1000
  
  dev.off()
  
  # browser()
  # write_rds(model_data, 'data/model_data_art.rds')
  
  return(
    cbind(
      suppression = threshold,
      slope_mean_link_identity = mean(model_data$slope_link_identity),
      pvalue_link_identity_slope_diff_zero = t.test((model_data$slope_link_identity))$p.value
    )
  )
}

ODn_regression_function_supressed <- function(dat, i=i) {
  # if(i==13)browser()
  model1 <- glm(
    formula = sedia_ODn_with_noise ~ time_on_ART + ART_status,
    family = gaussian(link = "identity"),
    data = dat
  )

  return(cbind(
    id = unique(dat$subject_label_blinded),
    intercept_link_identity = model1$coefficients[[1]],
    slope_link_identity = model1$coefficients[[2]]
  ))
}

results <- c()
for (i in c(1000)) {# 100, 400, 
  x <- slopes_for_suppressed(ODn_vl_data = full_dataset %>% 
                               # filter(time_on_ART <=365.25) %>%
                               mutate(ART_status = ifelse(time_on_ART > 0, 1, 2)), 
                             threshold = i)
  results <- rbind(results, x)
}


################################################################################
#' Describe Sedia LAg value slopes for records with last vl<1000 copies/ml to
#' the first visit where vl >1000 ('to peak') and vice versa (to tough)
################################################################################
data_intermitent_suppression_selected_visits <- read_delim("data/full_africos_data_with_ODn.csv", 
                                                         delim = ";", escape_double = FALSE, trim_ws = TRUE)  %>%
  filter(!is.na(to_peak)) %>%
  select(subject_label_blinded = subjid, #days_since_eddi, 
         Sex, Age, 
         test_date = Date_Collected, sedia_ODn = ODn, 
         viral_load = `Viral Load`, visits = Visit, to_peak) %>%
  left_join(
    readxl::read_excel("data/AFRICOS_Sempa Data Pull_24Jul23.xlsx", 
                       sheet = "Sheet1") %>%
      mutate(id = paste(Date_Collected, subjid, sep = '_')) %>%
      group_by(subjid) %>%
      mutate(switch_date = as.Date(as.numeric(ifelse(`ARV Start Date` == '.', '', `ARV Start Date`)), 
                                   origin = "1900-01-01"),
             art_start = as.Date(min(as.numeric(ifelse(`ARV Start Date` == '.', '', `ARV Start Date`)), na.rm = T), 
                                 origin = "1900-01-01")) %>%
      ungroup() %>%
      dplyr::select(subject_label_blinded = subjid, art_start) %>%
      distinct(subject_label_blinded, .keep_all = T), by = 'subject_label_blinded') %>%
  mutate(time_on_ART = as.numeric(as.Date(test_date, origin = "1900-01-01") - art_start),
         cohort = 'africos',
         # subject_label_blinded = as.double(cur_group_id()),
         Age = as.numeric(Age),
         viral_load = as.numeric(viral_load)) %>%
  select(subject_label_blinded, #days_since_eddi, 
         time_on_ART, Sex, Age, 
         test_date, sedia_ODn, 
         viral_load, visits, to_peak) %>%
  filter(time_on_ART >=0)


slopes_for_unsuppressed <- function(ODn_vl_data, threshold) {
  browser()
  dir.create(paste("unsupp", threshold, Sys.Date(), sep = "_"))
  folder_name <- paste("unsupp", threshold, Sys.Date(), sep = "_")
  
  # jpeg(paste(folder_name, "/", "Figure1_art", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)
  # 
  # print(
  #   ggplot(
  #     data = ODn_vl_data %>%
  #       # filter(suprressed_throughout_followup == 1) %>%
  #       group_by(subject_label_blinded),
  #     aes(x = time_on_ART, y = sedia_ODn)
  #   ) + # as.factor(id)
  #     geom_line(aes(color = as.factor(id)), size = 1) + # aes(color = id, linetype = Group), group = Group
  #     theme(
  #       text = element_text(size = 20),
  #       plot.title = element_text(hjust = 0.5),
  #       axis.line = element_line(colour = "black"),
  #       axis.text = element_text(size = 18),
  #       axis.title = element_text(size = 18),
  #       panel.background = element_blank(),
  #       panel.border = element_blank(),
  #       plot.margin = unit(c(0, 0, 0, 0), "null"),
  #       legend.position = "none"
  #     ) +
  #     ylab("Sedia LAg ODn") +
  #     xlab("Time on ART (Days)")
  # )
  # dev.off()
  
  model_data <- data.frame(
    id = NA, interceptlink_identity = NA, slope_link_identity = NA,
    intercept_link_log = NA, slope_link_log = NA
  )
  # set.seed(11)
  dat <- ODn_vl_data %>%
    filter(!is.na(time_on_ART)) %>%
    mutate(id = subject_label_blinded)
  
  counter <- 0
  for (i in 1:length(unique(dat$id))) {
    counter <- counter + 1
    model_data[counter, ] <- ODn_regression_function_unsupressed(
      dat = subset(dat, dat$id == unique(dat$id)[i])
    )
  }
  # browser()
  model_data <- model_data %>%
    mutate(`slope link identity` = as.numeric(slope_link_identity))
  
  hist_plot1 <- hist(model_data$`slope link identity`, breaks = 30, freq = FALSE, ylab = '', xlab = 'Individual slopes', main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  area_hist1 <- sum(hist_plot1$counts* abs(hist_plot1$mids[1]-hist_plot1$mids[2]))
  n_values1 <- sum(hist_plot1$counts / area_hist1)
  
  set.seed(11)
  x <- rnorm(n_values1, mean = mean(model_data$`slope link identity`, na.rm = T), sd = sd(model_data$`slope link identity`))
  dx <- density(x)

  hist_plot2 <- hist(model_data$`slope link log`, breaks = 30, freq = FALSE, ylab = '', xlab = 'Individual slopes', main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  area_hist2 <- sum(hist_plot2$counts* abs(hist_plot2$mids[1]-hist_plot2$mids[2]))
  n_values2 <- sum(hist_plot2$counts / area_hist2)
  set.seed(11)
  x1 <- rnorm(n_values2, mean = mean(model_data$`slope link log`, na.rm = T), sd = sd(model_data$`slope link log`))
  dx1 <- density(x1)
  jpeg(paste(folder_name, "/", "Figure_link_combined_art", threshold, ".jpeg", sep = ""), units = "in", width = 15, height = 8, res = 300)
    par(mfrow = c(1, 2))
    hist(model_data$`slope link identity`, breaks = 30, freq = FALSE, xlim = c(-0.02, 0.04), ylab = '', xlab = 'Slopes of individuals', main = 'A', col = 'black', cex.lab = 2, cex.axis = 2, xaxt="n", cex.main = 2) # , main = paste("Density plot link=identity_", threshold)
    axis(1, at = seq(-0.02, 0.04, by = 0.02), labels = seq(-0.02, 0.04, by = 0.02), cex.lab = 2, cex.axis = 2)
    lines(dx, lwd = 3, col = "red")
    
    hist(model_data$`slope link log`, breaks = 30, freq = FALSE, xlim = c(-0.020, 0.040), ylab = '', xlab = 'Slopes of individuals', main = 'B', col = 'black', cex.lab = 2, cex.axis = 2, xaxt="n", cex.main = 2) # , main = paste("Density plot link=identity_", threshold)
    axis(1, at = seq(-0.020, 0.040, by = 0.010), labels = seq(-0.020, 0.040, by = 0.010), cex.lab = 2, cex.axis = 2)
    lines(dx, lwd = 3, col = "red")
    dev.off()
write_rds(model_data, 'data/model_data_from_suppressed_art.rds')
  return(
    cbind(
      suppression = threshold,
      slope_mean_link_identity = mean(as.numeric(model_data$slope_link_identity)),
      pvalue_link_identity_slope_diff_zero = t.test((as.numeric(model_data$slope_link_identity)))$p.value
    )
  )
}

ODn_regression_function_unsupressed <- function(dat) {
  # browser()
  model1 <- glm(
    formula = sedia_ODn ~ time_on_ART,
    family = gaussian(link = "identity"),
    data = dat
  )
  
  return(cbind(
    id = unique(dat$id),
    intercept_link_identity = model1$coefficients[[1]],
    slope_link_identity = model1$coefficients[[2]]#,
    # intercept_link_log = model2$coefficients[[1]],
    # slope_link_log = model2$coefficients[[2]]
  ))
}

x_to_peak <- slopes_for_unsuppressed(
  ODn_vl_data = data_intermitent_suppression_selected_visits, #%>%
  #   filter(to_peak == 1) %>%
  #   unite("id", c(subject_label_blinded, set_to_peak)),
  threshold = 'to_unsuppressed'
)

results_unsuppressed <- rbind(x_to_peak, x_to_tough)

write.csv(results_unsuppressed, "output_table/results_unsuppressed.csv")

#########################################################################################
#'Using pooled standard deveiation
#'#######################################################################################

sigma_ODn_func <- function(data_set, threshold){
  # browser()
  data_generated <- data_set %>%
    select(
      subject_label_blinded, time_on_ART, test_date, sedia_ODn, viral_load # , art_initiation_date, aids_diagnosis_date,
      # art_interruption_date, art_resumption_date, treatment_naive,
      # on_treatment, first_treatment
    ) %>%
    arrange(subject_label_blinded, test_date) %>%
    # distinct(subject_label_blinded, test_date, .keep_all = T) %>%
    filter(!is.na(viral_load)) %>%
    group_by(subject_label_blinded) %>%
    mutate(flagvl = ifelse(viral_load <= threshold, 0, 1)) %>%
    mutate(suprressed_throughout_followup = ifelse(mean(flagvl) == 0, 1, 0)) %>%
    mutate(visits = 1:length(subject_label_blinded)) %>%
    mutate(n_visits = max(visits)) %>%
    filter(n_visits > 1) %>% # because at this stage we want to compare a with the group average
    ungroup() %>%
    filter(suprressed_throughout_followup == 1) %>%
    mutate(id = as.character(subject_label_blinded)) %>%
    dplyr::select(id, subject_label_blinded, sedia_ODn)%>%
    group_by(id) %>%
    summarise(n_length = n(),
              sd_value = sd(sedia_ODn))%>%
    mutate(numerator = ((n_length-1) * sd_value^2))
  g_stddev <- sum(data_generated$numerator)
  n_samples <- sum(data_generated$n_length)
  sigma_ODn <- (g_stddev / (n_samples-length(data_generated$n_length)))^.5
  
  return(sigma_ODn)
}
# sigma_ODn_func(sedia_generic, threshold = 1000)

compare_value_with_others <- function(data_set, threshold, sigma_ODn) {
  # browser()
  data_generated <- data_set %>%
    select(
      subject_label_blinded, time_on_ART, test_date, sedia_ODn, viral_load # , art_initiation_date, aids_diagnosis_date,
      # art_interruption_date, art_resumption_date, treatment_naive,
      # on_treatment, first_treatment
    ) %>%
    arrange(subject_label_blinded, test_date) %>%
    # distinct(subject_label_blinded, test_date, .keep_all = T) %>%
    filter(!is.na(viral_load)) %>%
    group_by(subject_label_blinded) %>%
    mutate(flagvl = ifelse(viral_load <= threshold, 0, 1)) %>%
    mutate(suprressed_throughout_followup = ifelse(mean(flagvl) == 0, 1, 0)) %>%
    mutate(visits = 1:length(subject_label_blinded)) %>%
    mutate(n_visits = max(visits)) %>%
    filter(n_visits > 1) %>% # because at this stage we want to compare a with the group average
    ungroup() %>%
    filter(suprressed_throughout_followup == 1) %>%
    mutate(id = as.character(subject_label_blinded)) %>%
    dplyr::select(id, subject_label_blinded, sedia_ODn, time_on_ART)
  
  results1 <- c() 
  ids <- unique(data_generated$id)
  # browser()
  for (i in 1:length(ids)) {
    ODn_values1 <- (data_generated %>%
                      filter(id == ids[i]) # %>%number_ids[i]
    )$sedia_ODn
    results_by_id <- c()
    for (j in 1:length(ODn_values1)) {
      test_value <- ODn_values1[j]
      list_values <- ODn_values1[!(ODn_values1 %in% test_value)]
      # z_test <- p_Value_z_test(list_values = list_values, test_value = test_value)
      length_list_value <- length(list_values)
      sigma_mu_s <- sigma_ODn / length_list_value^.5 #sd(list_values)
      mu_s <- mean(list_values)
      sigma_Y <- (sigma_mu_s^2  +  sigma_ODn^2)^.5#((((length(list_values) - 1) * sd_primary_set^2) + g_stddev) / ((n_samples + length_list_value) - (length(ids) - 1)))^.5
      z_test <- (test_value - mu_s) / sigma_Y
      set.seed(11)
      Z <- Normal(0, 1)  # make a standard normal r.v.
      p_value <- 1 - cdf(Z, abs(z_test)) + cdf(Z, -abs(z_test))
      results_by_id <- rbind(results_by_id, cbind(id = ids[i],value = test_value, z_stat = z_test, p_value = p_value)) # value = test_value,
    }
    results1 <- rbind(results1, results_by_id)
  }
  
  return(cbind(results1, threshold = threshold))
}

results_suppressed <- c()
for (i in c(1000)) {#100, 400, 
  # results1 <- c()
  results_suppressed <- rbind(
    results_suppressed, 
    compare_value_with_others(
      data_set = full_dataset %>% 
        filter(Group == 'early suppressions') %>%
        filter(time_on_ART >0 & viral_load < 1000), threshold = i,
      sigma_ODn = sigma_ODn_func(full_dataset %>% 
                                   filter(Group == 'early suppressions') %>%
                                   filter(time_on_ART >0 & viral_load < 1000), threshold = 1000)
      )
    )
  
}

##14660337, 22432139, 69689704,83908543, 99138302

results_suppressed1 <- as_tibble(results_suppressed) %>%
  distinct(id,z_stat, p_value, .keep_all = T) %>%
  # filter(threshold==1000) %>%
  mutate(significance = ifelse(as.numeric(p_value) < 0.05, TRUE, FALSE),
         `P-value` = as.numeric(p_value))

table(results_suppressed1$significance)
# hist(as.numeric(results_suppressed1$p_value), breaks = 50)

graph_results_suppressed1<-ggplot(results_suppressed1, aes(x=`P-value`)) +
  geom_histogram(color="black", fill="red", position="dodge")+
  geom_vline( aes(xintercept=0.05),
              linetype="dashed", size = 2)+
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    legend.position = "none"
  )
# graph_results_suppressed1
jpeg('other_figures/compare_value_with_others_art.jpeg', units = "in", width = 8, height = 6, res = 300)
graph_results_suppressed1
dev.off()

#' Compare the first LAg visit with VL > 1000 copies, with averages (backward moving averages) of the previous readings 
compare_lastvalue_with_previous <- function(data_set, sigma_ODn) {
  # browser()
  data_generated <- data_set %>% # data_intermitent_suppression_selected_visits
    filter(all_visits_to_peak == 1) %>%
    mutate(id = as.character(subject_label_blinded)) %>%
    dplyr::select(id, time_on_ART, test_date, sedia_ODn, viral_load, to_peak, all_visits_to_peak)
  results1 <- c()
  ids <- unique(data_generated$id)
  for (i in 1:length(ids)) {
    ODn_values1 <- (data_generated %>%
                      filter(id == ids[i])
    )$sedia_ODn
    results_by_id <- c()
    last_value_in_vec <- ODn_values1[length(ODn_values1)] # g <- 6
    new_vector <- ODn_values1[!(ODn_values1 %in% last_value_in_vec)]
    # browser()
    for (j in 1:(length(new_vector) - 1)) {
      x <- length(new_vector)
      y <- length(new_vector) - j
      list_values <- new_vector[x:y]
      test_value <- last_value_in_vec # ODn_values1[j]
      length_list_value <- length(list_values)
      sigma_mu_s <- sigma_ODn / length_list_value^.5 #sd(list_values)
      mu_s <- mean(list_values)
      # browser()
      sigma_Y <- (sigma_mu_s^2  +  sigma_ODn^2)^.5#((((length(list_values) - 1) * sd_primary_set^2) + g_stddev) / ((n_samples + length_list_value) - (length(ids) - 1)))^.5
      z_test <- (test_value - mu_s) / sigma_Y
      set.seed(11)
      Z <- Normal(0, 1)  # make a standard normal r.v.
      p_value <- 1 - cdf(Z, abs(z_test)) + cdf(Z, -abs(z_test))
      results_by_id <- rbind(results_by_id, cbind(id = ids[i],value = test_value, z_stat = z_test, p_value = p_value)) # value = test_value,
      
    }
    # print(ids[i])
    results1 <- rbind(results1, results_by_id)
  }
  return(results1)
}

compare_first_peak_value <- as_tibble(
  compare_lastvalue_with_previous(
    data_set = data_intermitent_suppression_selected_visits,
    sigma_ODn = sigma_ODn_func(sedia_generic, threshold = 1000))) %>%
  mutate(`significance 99%` = ifelse(as.numeric(z_stat) > qnorm(0.99), TRUE, FALSE),
         `significance 98.5%` = ifelse(as.numeric(z_stat) > qnorm(0.985), TRUE, FALSE),
         `significance 98%` = ifelse(as.numeric(z_stat) > qnorm(0.98), TRUE, FALSE),
         `significance 97.5%` = ifelse(as.numeric(z_stat) > qnorm(0.975), TRUE, FALSE),
         `significance 95%` = ifelse(as.numeric(z_stat) > qnorm(0.95), TRUE, FALSE),
         `P-value` = as.numeric(p_value))
