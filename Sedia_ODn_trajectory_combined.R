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
final_dataset <- read_csv("Sempa_final_pull_with_results.csv") %>%
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
  filter(Group == 'early suppressions') %>%
  select(subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load, visits, baseline_visit) 

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

sd_sedia_ODn_vs_ODn_bc

noise <- summary(glm(sd_Sedia_ODn ~ mean_Sedia_ODn, data = sedia_distribution_blinded_controls %>%
  mutate(sd_Sedia_ODn = `sigma Sedia ODn`, mean_Sedia_ODn = `mean Sedia ODn`)))

sedia_generic <- sedia_generic %>%
  filter(visit_hivstatus == "P") %>%
  # rename(sedia_ODn = result...72) %>%
  filter(visit_id != 21773 & visit_id != 21783 & visit_id != 21785)

######################################################
#' Analysis among VL suppressed at 100, 400, and 1000
######################################################

slopes_for_suppressed <- function(ODn_vl_data, threshold) {
  # browser()
  dir.create(paste("supp_", threshold, Sys.Date(), sep = "_"))
  folder_name <- paste("supp_", threshold, Sys.Date(), sep = "_")
  data_generated <- ODn_vl_data %>%
    select(
      subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load # , art_initiation_date, aids_diagnosis_date,
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
      aes(x = days_since_eddi, y = sedia_ODn)
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
      xlab("Time since infection (Days)") +
      ggtitle(paste("Sedia ODn vs EDDI_", threshold))
  )
  dev.off()
# browser()
  model_data <- data.frame(
    id = NA, interceptlink_identity = NA, slope_link_identity = NA,
    intercept_link_log = NA, slope_link_log = NA
  )
  set.seed(11)
  dat <- data_generated %>%
    filter(suprressed_throughout_followup == 1) %>%
    mutate(
      sedia_ODn_with_noise = sedia_ODn + rnorm(n = length(subject_label_blinded), mean = 0, sd = noise$coefficients[1, 1] + noise$coefficients[2, 1] * sedia_ODn),
      sedia_ODn_with_noise = ifelse(sedia_ODn_with_noise < 0, 0, sedia_ODn_with_noise),
      sedia_ODn_with_noise = ifelse(is.na(sedia_ODn_with_noise), sedia_ODn, sedia_ODn_with_noise)
    ) %>%
    filter(!is.na(days_since_eddi)) %>%
    select(subject_label_blinded, test_date, days_since_eddi, sedia_ODn, sedia_ODn_with_noise, viral_load)
  counter <- 0
  for (i in 1:length(unique(dat$subject_label_blinded))) {
    counter <- counter + 1
    model_data[counter, ] <- ODn_regression_function_supressed(
      dat = subset(dat, dat$subject_label_blinded == unique(dat$subject_label_blinded)[i])
    )
    # print(i)
  }
  # summary(model_data$slope_link_identity, na.rm = T)
  # browser()
  hist_plot1 <- hist(model_data$slope_link_identity, breaks = 30, freq = FALSE, ylab = '', xlab = 'Individual slopes', main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  area_hist1 <- sum(hist_plot1$counts* abs(hist_plot1$mids[1]-hist_plot1$mids[2]))
  n_values1 <- sum(hist_plot1$counts / area_hist1)
  set.seed(11)
  x <- rnorm(n_values1, mean = mean(model_data$slope_link_identity, na.rm = T), sd = sd(model_data$slope_link_identity))
  dx <- density(x)
  # jpeg(paste(folder_name, "/", "Figure_link_identity", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)
  # hist(model_data$slope_link_identity, breaks = 30, freq = FALSE, ylab = '', xlab = 'Slopes of individuals', main = '', col = 'black', cex.lab = 2, cex.axis = 2) #, main = paste("Density plot link=identity_", threshold)
  # 
  # # hist(model_data$slope_link_identity, breaks = 30, freq = FALSE, ylab = '', xlab = 'Individual slopes', main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  # lines(dx, lwd = 3, col = "red")
  # 
  # dev.off()
  # summary(model_data$slope_link_log, na.rm = T)
  hist_plot2 <- hist(model_data$slope_link_log, breaks = 30, freq = FALSE, ylab = '', xlab = 'Individual slopes', main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  area_hist2 <- sum(hist_plot2$counts* abs(hist_plot2$mids[1]-hist_plot2$mids[2]))
  n_values2 <- sum(hist_plot2$counts / area_hist2)
  set.seed(11)
  x1 <- rnorm(n_values2, mean = mean(model_data$slope_link_log, na.rm = T), sd = sd(model_data$slope_link_log))
  dx1 <- density(x1)
  
  # jpeg(paste(folder_name, "/", "Figure_link_log", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)
  # 
  # hist(model_data$slope_link_log, breaks = 30, freq = FALSE, ylab = '', xlab = 'Slopes of individuals', main = '', col = 'black', cex.lab = 2, cex.axis = 2) # , main = paste("Density plot link=identity_", threshold)
  # lines(dx, lwd = 3, col = "red")
  # 
  # dev.off()
  # browser()#
  jpeg(paste(folder_name, "/", "Figure_combined", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 8, res = 300)
  ggplot(
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
  # hist(model_data$slope_link_identity, breaks = 30, xlim = c(-0.010, 0.010), xaxt="n", ylab = 'Counts', xlab = 'Slopes of individuals', main = '', col = 'black', cex.lab = 2, cex.axis = 2) #, main = paste("Density plot link=identity_", threshold)
  # axis(1, at = seq(-.010, .01, by = .02), labels = seq(-.010, .01,by = .02) )
  # hist(model_data$slope_link_log, breaks = 30, ylab = 'Counts', xlab = 'Slopes of individuals', main = 'B', col = 'black', cex.lab = 2, cex.axis = 2) # , main = paste("Density plot link=identity_", threshold)
  # lines(dx1, lwd = 3, col = "red")
  dev.off()
  
  # browser()
  write_rds(model_data, 'data/model_data.rds')
  
  return(
    cbind(
      suppression = threshold,
      slope_mean_link_identity = mean(model_data$slope_link_identity),
      pvalue_link_identity_slope_diff_zero = t.test((model_data$slope_link_identity))$p.value,
      slope_mean_link_log = mean(model_data$slope_link_log),
      pvalue_link_log_slope_diff_zero = t.test((model_data$slope_link_log))$p.value
    )
  )
}

ODn_regression_function_supressed <- function(dat) {
  # browser()
  model1 <- glm(
    formula = sedia_ODn_with_noise ~ days_since_eddi,
    family = gaussian(link = "identity"),
    data = dat
  )

  model2 <- glm(
    formula = sedia_ODn_with_noise ~ days_since_eddi,
    family = gaussian(link = "log"),
    data = dat
  )

  return(cbind(
    id = unique(dat$subject_label_blinded),
    intercept_link_identity = model1$coefficients[[1]],
    slope_link_identity = model1$coefficients[[2]],
    intercept_link_log = model2$coefficients[[1]],
    slope_link_log = model2$coefficients[[2]]
  ))
}


results <- c()
for (i in c(1000)) {# 100, 400, 
  x <- slopes_for_suppressed(ODn_vl_data = sedia_generic, threshold = i)
  results <- rbind(results, x)
}
write.csv(results, "output_table/results_suppressed.csv") # results
set.seed(11);Z <- Normal(0, 1)
model_data <- read_rds('data/model_data.rds') %>%
  mutate(sigma_slope_identity = (sd(slope_link_identity)/(length(slope_link_identity) - 1))^0.5,
         sigma_slope_log = (sd(slope_link_log)/(length(slope_link_log) - 1))^0.5) %>%
  group_by(id) %>%
  mutate(z_value_identity = (slope_link_identity - 0) / sigma_slope_identity,
         `P values identity` = 1 - cdf(Z, abs(z_value_identity)) + cdf(Z, -abs(z_value_identity)),
         z_value_log = (slope_link_log - 0) / sigma_slope_log,
         `P values log` = 1 - cdf(Z, abs(z_value_log)) + cdf(Z, -abs(z_value_log))
         )

################################################################
#' Analysis among early VL suppressed test straddling ART start
################################################################

slopes_for_suppressed_straddling_ARTstart <- function(ODn_vl_data, threshold) {
  # browser()
  dir.create(paste("strad_ARTstart_supp_", threshold, Sys.Date(), sep = "_"))
  folder_name <- paste("strad_ARTstart_supp_", threshold, Sys.Date(), sep = "_")
  data_generated <- ODn_vl_data# %>%
    # select(
    #   subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load # , art_initiation_date, aids_diagnosis_date,
    #   # art_interruption_date, art_resumption_date, treatment_naive,
    #   # on_treatment, first_treatment
    # ) %>%
    # arrange(subject_label_blinded, test_date) %>%
    # # distinct(subject_label_blinded, test_date, .keep_all = T) %>%
    # filter(!is.na(viral_load)) %>%
    # group_by(subject_label_blinded) %>%
    # mutate(flagvl = ifelse(viral_load <= threshold, 0, 1)) %>%
    # mutate(suprressed_throughout_followup = ifelse(mean(flagvl) == 0, 1, 0)) %>%
    # mutate(visits = 1:length(subject_label_blinded)) %>%
    # mutate(n_visits = max(visits)) %>%
    # filter(n_visits > 1) %>%
    # ungroup()
  
  jpeg(paste(folder_name, "/", "Figure1", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)
  
  print(
    ggplot(
      data = data_generated %>%
        # filter(suprressed_throughout_followup == 1) %>%
        group_by(subject_label_blinded),
      aes(x = days_since_eddi, y = sedia_ODn)
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
      xlab("Time since infection (Days)") +
      ggtitle(paste("Sedia ODn vs EDDI_", threshold))
  )
  dev.off()
  
  model_data <- data.frame(
    id = NA, interceptlink_identity = NA, slope_link_identity = NA,
    intercept_link_log = NA, slope_link_log = NA
  )
  # browser()
  set.seed(11)
  dat <- data_generated %>%
    # filter(suprressed_throughout_followup == 1) %>%
    mutate(
      sedia_ODn_with_noise = sedia_ODn + rnorm(n = length(subject_label_blinded), mean = 0, sd = noise$coefficients[1, 1] + noise$coefficients[2, 1] * sedia_ODn),
      sedia_ODn_with_noise = ifelse(sedia_ODn_with_noise < 0, 0, ifelse(viral_load>999, sedia_ODn, sedia_ODn_with_noise)),
      sedia_ODn_with_noise = ifelse(is.na(sedia_ODn_with_noise), sedia_ODn, sedia_ODn_with_noise)
    ) %>%
    filter(!is.na(days_since_eddi)) %>%
    select(subject_label_blinded, test_date, days_since_eddi, sedia_ODn, sedia_ODn_with_noise, viral_load)
  # browser()
  counter <- 0
  for (i in 1:length(unique(dat$subject_label_blinded))) {
    counter <- counter + 1
    # print(unique(dat$subject_label_blinded)[i])
    model_data[counter, ] <- ODn_regression_function_supressed_straddling_ARTstart(
      dat = subset(dat, dat$subject_label_blinded == unique(dat$subject_label_blinded)[i])
    )
    # print(i)
  }
  # summary(model_data$slope_link_identity, na.rm = T)
  # browser()
  hist_plot1 <- hist(model_data$slope_link_identity, breaks = 30, freq = FALSE, ylab = '', xlab = 'Individual slopes', main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  area_hist1 <- sum(hist_plot1$counts* abs(hist_plot1$mids[1]-hist_plot1$mids[2]))
  n_values1 <- sum(hist_plot1$counts / area_hist1)
  set.seed(11)
  x <- rnorm(n_values1, mean = mean(model_data$slope_link_identity, na.rm = T), sd = sd(model_data$slope_link_identity))
  dx <- density(x)
  # jpeg(paste(folder_name, "/", "Figure_link_identity", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)
  # hist(model_data$slope_link_identity, breaks = 30, freq = FALSE, ylab = '', xlab = 'Slopes of individuals', main = '', col = 'black', cex.lab = 2, cex.axis = 2) # , main = paste("Density plot link=identity_", threshold)
  # 
  # # hist(model_data$slope_link_identity, breaks = 30, freq = FALSE, main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  # lines(dx, lwd = 3, col = "red")
  # dev.off()
  # summary(model_data$slope_link_log, na.rm = T)
  hist_plot2 <- hist(model_data$slope_link_log, breaks = 30, freq = FALSE, ylab = '', xlab = 'Individual slopes', main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  area_hist2 <- sum(hist_plot2$counts* abs(hist_plot2$mids[1]-hist_plot2$mids[2]))
  n_values2 <- sum(hist_plot2$counts / area_hist2)
  set.seed(11)
  x1 <- rnorm(n_values2, mean = mean(model_data$slope_link_log, na.rm = T), sd = sd(model_data$slope_link_log))
  dx1 <- density(x1)
  # jpeg(paste(folder_name, "/", "Figure_link_log", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)
  # hist(model_data$slope_link_log, breaks = 30, freq = FALSE, ylab = '', xlab = 'Slopes of individuals', main = '', col = 'black', cex.lab = 2, cex.axis = 2) # main = paste("Density plot link=identity_", threshold), 
  # # hist(model_data$slope_link_log, breaks = 30, freq = FALSE, main = paste("Mean slope link=log_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  # lines(dx, lwd = 3, col = "red")
  # 
  # print(
  #   model_data %>%
  #     ggplot(aes(x = slope_link_log)) +
  #     geom_histogram(color = "#e9ecef", alpha = 0.6, bins = 30, position = "identity") +
  #     scale_fill_manual(values = c("#69b3a2", "#404080")) +
  #     labs(fill = "") +
  #     geom_density(
  #       lwd = 1.2,
  #       linetype = 1,
  #       colour = 2
  #     ) +
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
  #     ggtitle(paste("histogram mean slope link=log_", threshold))
  # )
  # dev.off()
  
  jpeg(paste(folder_name, "/", "Figure_link_combined", threshold, ".jpeg", sep = ""), units = "in", width = 15, height = 8, res = 300)
  par(mfrow = c(1, 2))
  hist(model_data$slope_link_identity, breaks = 30, freq = FALSE, ylab = '', xlab = 'Slopes of individuals', main = 'A', col = 'black', cex.lab = 2, cex.axis = 2) # , main = paste("Density plot link=identity_", threshold)
  lines(dx, lwd = 3, col = "red")
  hist(model_data$slope_link_log, breaks = 30, freq = FALSE, ylab = '', xlab = 'Slopes of individuals', main = 'B', col = 'black', cex.lab = 2, cex.axis = 2) # main = paste("Density plot link=identity_", threshold), 
  lines(dx, lwd = 3, col = "red")
  dev.off()
  
  
  return(
    cbind(
      # suppression = threshold,
      slope_mean_link_identity = mean(model_data$slope_link_identity),
      pvalue_link_identity_slope_diff_zero = t.test((model_data$slope_link_identity))$p.value,
      slope_mean_link_log = mean(model_data$slope_link_log),
      pvalue_link_log_slope_diff_zero = t.test((model_data$slope_link_log))$p.value
    )
  )
}

ODn_regression_function_supressed_straddling_ARTstart <- function(dat) {
  # browser()
  model1 <- glm(
    formula = sedia_ODn_with_noise ~ days_since_eddi,
    family = gaussian(link = "identity"),
    data = dat
  )
  
  model2 <- glm(
    formula = sedia_ODn_with_noise ~ days_since_eddi,
    family = gaussian(link = "log"),
    data = dat,
    start = c(model1$coefficients[[1]], model1$coefficients[[2]])
  )
  
  return(cbind(
    id = unique(dat$subject_label_blinded),
    intercept_link_identity = model1$coefficients[[1]],
    slope_link_identity = model1$coefficients[[2]],
    intercept_link_log = model2$coefficients[[1]],
    slope_link_log = model2$coefficients[[2]]
  ))
}


results <- c()
# for (i in c(100, 400, 1000)) {
  x <- slopes_for_suppressed_straddling_ARTstart(ODn_vl_data = final_dataset, threshold = 'ES')
  results <- rbind(results, x)
# }
write.csv(results, "output_table/results_suppressed_straddling_ARTstart.csv") # results


################################################################################
#' Describe Sedia LAg value slopes for records with last vl<1000 copies/ml to
#' the first visit where vl >1000 ('to peak') and vice versa (to tough)
################################################################################
data_intermitent_suppression <- read_csv("data/20180410-EP-LAgSedia-Generic.csv") %>%
  mutate(sedia_ODn = `result...15`) %>%
  select(
    subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load # , art_initiation_date, aids_diagnosis_date,
    # art_interruption_date, art_resumption_date, treatment_naive,
    # on_treatment, first_treatment
  ) %>%
  arrange(subject_label_blinded, test_date) %>%
  # distinct(subject_label_blinded, test_date, .keep_all = T) %>%
  filter(!is.na(viral_load)) %>%
  group_by(subject_label_blinded) %>%
  mutate(flagvl = ifelse(viral_load > 1000, 0, 1)) %>%
  mutate(irregular_suprression = ifelse(mean(flagvl) == 0, 1, 0)) %>%
  filter(irregular_suprression == 0) %>%
  mutate(visits = 1:length(subject_label_blinded)) %>%
  mutate(n_visits = max(visits)) %>%
  filter(n_visits > 1) %>%
  arrange(subject_label_blinded, days_since_eddi) %>%
  ungroup()


#' does the slope for ODn among those with ODn<2 differ from those with ODn>2 in
#' patients peaking?
# write.csv(data_intermitent_suppression, "output_table/intermitent_suppression.csv")
data_intermitent_suppression_selected_visits <- read.csv("output_table/intermitent_suppression_selected.csv") %>%
  mutate(Cohort = 'CEPHIA') %>%
  bind_rows(read.csv('data/africos_data_with_ODn.csv') %>%
              mutate(Cohort = 'AFRICOS')) %>%
  mutate(all_visits_to_peak = ifelse(subject_label_blinded == 28848447, NA, all_visits_to_peak))
data_intermitent_suppression_toPeak_visits <- data_intermitent_suppression_selected_visits %>%
  filter(!is.na(to_peak)) %>%
  unite("id", c(subject_label_blinded, set_to_peak))

data_intermitent_suppression_totrough_visits <- data_intermitent_suppression_selected_visits %>%
  filter(!is.na(to_trough)) %>%
  unite("id", c(subject_label_blinded, set_to_trough))
summary(to_peak <- data_intermitent_suppression_selected_visits %>%
  filter(!is.na(to_peak)) %>%
  unite("id", c(subject_label_blinded, set_to_peak)) %>%
  group_by(id) %>%
  summarise(x = mean(diff(days_since_eddi))))

summary(to_trough <- data_intermitent_suppression_selected_visits %>%
  filter(!is.na(to_trough)) %>%
  unite("id", c(subject_label_blinded, set_to_trough)) %>%
  group_by(id) %>%
  summarise(x = mean(diff(days_since_eddi))))

slopes_for_unsuppressed <- function(ODn_vl_data, threshold) {
  # browser()
  dir.create(paste("unsupp", threshold, Sys.Date(), sep = "_"))
  folder_name <- paste("unsupp", threshold, Sys.Date(), sep = "_")

  jpeg(paste(folder_name, "/", "Figure1", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)

  print(
    ggplot(
      data = ODn_vl_data %>%
        # filter(suprressed_throughout_followup == 1) %>%
        group_by(id),
      aes(x = days_since_eddi, y = sedia_ODn)
    ) + # as.factor(id)
      geom_line(aes(color = as.factor(id)), size = 1) + # aes(color = id, linetype = Group), group = Group
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
      xlab("Time since infection (Days)")
  )
  dev.off()

  model_data <- data.frame(
    id = NA, interceptlink_identity = NA, slope_link_identity = NA,
    intercept_link_log = NA, slope_link_log = NA
  )
  # set.seed(11)
  dat <- ODn_vl_data %>%
    filter(!is.na(days_since_eddi))

  counter <- 0
  for (i in 1:length(unique(dat$id))) {
    counter <- counter + 1
    model_data[counter, ] <- ODn_regression_function_unsupressed(
      dat = subset(dat, dat$id == unique(dat$id)[i])
    )
  }
# browser()
  model_data <- model_data %>%
    mutate(`slope link identity` = as.numeric(slope_link_identity),
           `slope link log` = as.numeric(slope_link_log))
  
  hist_plot1 <- hist(model_data$`slope link identity`, breaks = 30, freq = FALSE, ylab = '', xlab = 'Individual slopes', main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  area_hist1 <- sum(hist_plot1$counts* abs(hist_plot1$mids[1]-hist_plot1$mids[2]))
  n_values1 <- sum(hist_plot1$counts / area_hist1)
  
  set.seed(11)
  x <- rnorm(n_values1, mean = mean(model_data$`slope link identity`, na.rm = T), sd = sd(model_data$`slope link identity`))
  dx <- density(x)
  # browser()
  # jpeg(paste(folder_name, "/", "Figure_link_identity", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)
  # 
  # hist(model_data$`slope link identity`, breaks = 30, freq = FALSE, xlim = c(-0.02, 0.04), ylab = '', xlab = 'Slopes of individuals', main = '', col = 'black', cex.lab = 2, cex.axis = 2, xaxt="n") # , main = paste("Density plot link=identity_", threshold)
  # axis(1, at = seq(-0.02, 0.04, by = 0.02), labels = seq(-0.02, 0.04, by = 0.02), cex.lab = 2, cex.axis = 2)
  # lines(dx, lwd = 3, col = "red")
  # dev.off()
  
  hist_plot2 <- hist(model_data$`slope link log`, breaks = 30, freq = FALSE, ylab = '', xlab = 'Individual slopes', main = paste("Mean slope link=identity_", threshold), col = 'black', cex.lab = 1.5, cex.axis = 1.5)
  area_hist2 <- sum(hist_plot2$counts* abs(hist_plot2$mids[1]-hist_plot2$mids[2]))
  n_values2 <- sum(hist_plot2$counts / area_hist2)
  set.seed(11)
  x1 <- rnorm(n_values2, mean = mean(model_data$`slope link log`, na.rm = T), sd = sd(model_data$`slope link log`))
  dx1 <- density(x1)
  
  
  # jpeg(paste(folder_name, "/", "Figure_link_log",threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)
  # 
  # hist(model_data$`slope link log`, breaks = 30, freq = FALSE, xlim = c(-0.020, 0.040), ylab = '', xlab = 'Slopes of individuals', main = '', col = 'black', cex.lab = 2, cex.axis = 2, xaxt="n") # , main = paste("Density plot link=identity_", threshold)
  # axis(1, at = seq(-0.020, 0.040, by = 0.010), labels = seq(-0.020, 0.040, by = 0.010), cex.lab = 2, cex.axis = 2)
  # lines(dx, lwd = 3, col = "red")
  # dev.off()
  
  #####
  if (threshold == 'to_unsuppressed'){
  jpeg(paste(folder_name, "/", "Figure_link_combined", threshold, ".jpeg", sep = ""), units = "in", width = 15, height = 8, res = 300)
  par(mfrow = c(1, 2))
  hist(model_data$`slope link identity`, breaks = 30, freq = FALSE, xlim = c(-0.02, 0.04), ylab = '', xlab = 'Slopes of individuals', main = 'A', col = 'black', cex.lab = 2, cex.axis = 2, xaxt="n", cex.main = 2) # , main = paste("Density plot link=identity_", threshold)
  axis(1, at = seq(-0.02, 0.04, by = 0.02), labels = seq(-0.02, 0.04, by = 0.02), cex.lab = 2, cex.axis = 2)
  lines(dx, lwd = 3, col = "red")
  
  hist(model_data$`slope link log`, breaks = 30, freq = FALSE, xlim = c(-0.020, 0.040), ylab = '', xlab = 'Slopes of individuals', main = 'B', col = 'black', cex.lab = 2, cex.axis = 2, xaxt="n", cex.main = 2) # , main = paste("Density plot link=identity_", threshold)
  axis(1, at = seq(-0.020, 0.040, by = 0.010), labels = seq(-0.020, 0.040, by = 0.010), cex.lab = 2, cex.axis = 2)
  lines(dx, lwd = 3, col = "red")
  dev.off()}
  
  if (threshold == 'from_unsuppressed'){
    jpeg(paste(folder_name, "/", "Figure_link_combined", threshold, ".jpeg", sep = ""), units = "in", width = 15, height = 8, res = 300)
    par(mfrow = c(1, 2))
    hist(model_data$`slope link identity`, breaks = 30, freq = FALSE, xlim = c(-0.03, 0.02), ylab = '', xlab = 'Slopes of individuals', main = 'A', col = 'black', cex.lab = 2, cex.axis = 2, xaxt="n", cex.main = 2) # , main = paste("Density plot link=identity_", threshold)
    axis(1, at = seq(-0.03, 0.02, by = 0.01), labels = seq(-0.03, 0.02, by = 0.01), cex.lab = 2, cex.axis = 2)
    lines(dx, lwd = 3, col = "red")
    
    hist(model_data$`slope link log`, breaks = 30, freq = FALSE, xlim = c(-0.040, 0.030), ylab = '', xlab = 'Slopes of individuals', main = 'B', col = 'black', cex.lab = 2, cex.axis = 2, xaxt="n", cex.main = 2) # , main = paste("Density plot link=identity_", threshold)
    axis(1, at = seq(-0.040, 0.030, by = 0.010), labels = seq(-0.040, 0.030, by = 0.010), cex.lab = 2, cex.axis = 2)
    lines(dx, lwd = 3, col = "red")
    dev.off()}
  
  if(threshold =='to_unsuppressed'){ write_rds(model_data, 'data/model_data_from_suppressed.rds')}
  write_rds(model_data, 'data/model_data_to_suppressed.rds')
  return(
    cbind(
      suppression = threshold,
      slope_mean_link_identity = mean(as.numeric(model_data$slope_link_identity)),
      pvalue_link_identity_slope_diff_zero = t.test((as.numeric(model_data$slope_link_identity)))$p.value,
      slope_mean_link_log = mean(as.numeric(model_data$slope_link_log)),
      pvalue_link_log_slope_diff_zero = t.test((as.numeric(model_data$slope_link_log)))$p.value
    )
  )
}

ODn_regression_function_unsupressed <- function(dat) {
  # browser()
  model1 <- glm(
    formula = sedia_ODn ~ days_since_eddi,
    family = gaussian(link = "identity"),
    data = dat
  )

  model2 <- glm(
    formula = sedia_ODn ~ days_since_eddi,
    family = gaussian(link = "log"),
    data = dat
  )

  return(cbind(
    id = unique(dat$id),
    intercept_link_identity = model1$coefficients[[1]],
    slope_link_identity = model1$coefficients[[2]],
    intercept_link_log = model2$coefficients[[1]],
    slope_link_log = model2$coefficients[[2]]
  ))
}

x_to_peak <- slopes_for_unsuppressed(
  ODn_vl_data = data_intermitent_suppression_selected_visits %>%
    filter(to_peak == 1) %>%
    unite("id", c(subject_label_blinded, set_to_peak)),
  threshold = 'to_unsuppressed'
)
x_to_tough <- slopes_for_unsuppressed(
  ODn_vl_data = data_intermitent_suppression_selected_visits %>%
    filter(to_trough == 1) %>%
    unite("id", c(subject_label_blinded, set_to_trough)),
  threshold = 'from_unsuppressed'
)
results_unsuppressed <- rbind(x_to_peak, x_to_tough)

write.csv(results_unsuppressed, "output_table/results_unsuppressed.csv")

#########################################################################################
#'Using pooled standard deveiation
#'#######################################################################################

sigma_ODn_func <- function(data_set, threshold){
  data_generated <- data_set %>%
    select(
      subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load # , art_initiation_date, aids_diagnosis_date,
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
  data_generated <- data_set %>%
    select(
      subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load # , art_initiation_date, aids_diagnosis_date,
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
    dplyr::select(id, subject_label_blinded, sedia_ODn)
  
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
for (i in c(100, 400,1000)) {# 
  # results1 <- c()
  results_suppressed <- rbind(results_suppressed, compare_value_with_others(data_set = sedia_generic, threshold = i,
                              sigma_ODn = sigma_ODn_func(sedia_generic, threshold = 1000)))
  
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
jpeg('other_figures/compare_value_with_others.jpeg', units = "in", width = 8, height = 6, res = 300)
graph_results_suppressed1
dev.off()
#' Compare the first LAg visit with VL > 1000 copies, with averages (backward moving averages) of the previous readings 
compare_lastvalue_with_previous <- function(data_set, sigma_ODn) {
  # browser()
  data_generated <- data_set %>% # data_intermitent_suppression_selected_visits
    filter(all_visits_to_peak == 1) %>%
    mutate(id = as.character(subject_label_blinded)) %>%
    dplyr::select(id, days_since_eddi, test_date, sedia_ODn, viral_load, to_peak, all_visits_to_peak)
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

compare_first_peak_value <- as_tibble(compare_lastvalue_with_previous(data_set = data_intermitent_suppression_selected_visits,
                                                                      sigma_ODn = sigma_ODn_func(sedia_generic, threshold = 1000))) %>%
  mutate(`significance 99%` = ifelse(as.numeric(z_stat) > qnorm(0.99), TRUE, FALSE),
         `significance 98.5%` = ifelse(as.numeric(z_stat) > qnorm(0.985), TRUE, FALSE),
         `significance 98%` = ifelse(as.numeric(z_stat) > qnorm(0.98), TRUE, FALSE),
         `significance 97.5%` = ifelse(as.numeric(z_stat) > qnorm(0.975), TRUE, FALSE),
         `significance 95%` = ifelse(as.numeric(z_stat) > qnorm(0.95), TRUE, FALSE),
         `P-value` = as.numeric(p_value))

# table(compare_first_peak_value$`significance 95%`)
# table(compare_first_peak_value$`significance 97.5%`)
# table(compare_first_peak_value$`significance 98%`)
# table(compare_first_peak_value$`significance 98.5%`)
# table(compare_first_peak_value$`significance 99%`)

sensitivity_value <- cbind(`Accuracy values` = c(
  round(table(compare_first_peak_value$`significance 95%`)[[2]] / (table(compare_first_peak_value$`significance 95%`)[[1]] + table(compare_first_peak_value$`significance 95%`)[[2]]), 3),
  round(table(compare_first_peak_value$`significance 97.5%`)[[2]] / (table(compare_first_peak_value$`significance 97.5%`)[[1]] + table(compare_first_peak_value$`significance 97.5%`)[[2]]), 3),
  round(table(compare_first_peak_value$`significance 98%`)[[2]] / (table(compare_first_peak_value$`significance 98%`)[[1]] + table(compare_first_peak_value$`significance 98%`)[[2]]), 3),
  round(table(compare_first_peak_value$`significance 98.5%`)[[2]] / (table(compare_first_peak_value$`significance 98.5%`)[[1]] + table(compare_first_peak_value$`significance 98.5%`)[[2]]), 3),
  round(table(compare_first_peak_value$`significance 99%`)[[2]] / (table(compare_first_peak_value$`significance 99%`)[[1]] + table(compare_first_peak_value$`significance 99%`)[[2]]), 3)
), level = c('Significance 95%', 'Significance 97.5%', 'Significance 98%', 'Significance 98.5%', 'Significance 99%'), Accuracy = rep('Sensitivity', 5)
)
graph_compare_first_peak_value<-ggplot(compare_first_peak_value, aes(x=`P-value`)) +
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
# graph_compare_first_peak_value
jpeg('other_figures/compare_lastvalue_with_previous.jpeg', units = "in", width = 8, height = 6, res = 300)
graph_compare_first_peak_value
dev.off()
#'comparing LAg value with the average of the previous LAg readings
compare_value_with_previous <- function(data_set, threshold, sigma_ODn) {
  # browser()
  data_generated <- data_set %>%
    select(
      subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load # , art_initiation_date, aids_diagnosis_date,
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
    dplyr::select(id, subject_label_blinded, sedia_ODn)
  ids <- unique(data_generated$id)
  results1 <- c()
  for (i in 1:length(ids)) {
    ODn_values1 <- (data_generated %>%
      filter(id == ids[i]) # %>%number_ids[i]
    )$sedia_ODn
    results_by_id <- c()
    for (j in 2:length(ODn_values1)) {
      test_value <- ODn_values1[j]
      list_values <- ODn_values1[(j - 1):1]
      length_list_value <- length(list_values)
      sigma_mu_s <- sigma_ODn / length_list_value^.5 #sd(list_values)
      mu_s <- mean(list_values)
      sigma_Y <- (sigma_mu_s^2  +  sigma_ODn^2)^.5#((((length(list_values) - 1) * sd_primary_set^2) + g_stddev) / ((n_samples + length_list_value) - (length(ids) - 1)))^.5
      z_test <- (test_value - mu_s) / sigma_Y
      set.seed(11)
      Z <- Normal(0, 1) # make a standard normal r.v.
      p_value <- 1 - cdf(Z, abs(z_test)) + cdf(Z, -abs(z_test))
      results_by_id <- rbind(results_by_id, cbind(id = ids[i], value = test_value, z_stat = z_test, p_value = p_value)) # value = test_value,
      
    }
    results1 <- rbind(results1, results_by_id)
  }

  return(cbind(results1, threshold = threshold))
}

results_suppressed <- c()
for (i in c(100, 400, 1000)) {
  # results1 <- c()
  results_suppressed <- rbind(results_suppressed, compare_value_with_previous(data_set = sedia_generic, threshold = i,
                                                                              sigma_ODn = sigma_ODn_func(sedia_generic, threshold = 1000)))
  
}
results_suppressed2 <- as_tibble(results_suppressed) %>%
  distinct(id, value, z_stat, p_value, .keep_all = T) %>%
  # filter(threshold==1000) %>%
  mutate(`significance 99%` = ifelse(as.numeric(z_stat) > qnorm(0.99), TRUE, FALSE),
         `significance 98.5%` = ifelse(as.numeric(z_stat) > qnorm(0.985), TRUE, FALSE),
         `significance 98%` = ifelse(as.numeric(z_stat) > qnorm(0.98), TRUE, FALSE),
         `significance 97.5%` = ifelse(as.numeric(z_stat) > qnorm(0.975), TRUE, FALSE),
         `significance 95%` = ifelse(as.numeric(z_stat) > qnorm(0.95), TRUE, FALSE),
         `P-value` = as.numeric(p_value)) %>%
  group_by(id) %>%
  mutate(n_comparison = 2:(length(id)+1))

#table(results_suppressed2$`significance 97.5%`)
specificity_value <- cbind(`Accuracy values` = c(
  round(table(results_suppressed2$`significance 95%`)[[1]] / (table(results_suppressed2$`significance 95%`)[[1]] + table(results_suppressed2$`significance 95%`)[[2]]), 3),
  round(table(results_suppressed2$`significance 97.5%`)[[1]] / (table(results_suppressed2$`significance 97.5%`)[[1]] + table(results_suppressed2$`significance 97.5%`)[[2]]), 3),
  round(table(results_suppressed2$`significance 98%`)[[1]] / (table(results_suppressed2$`significance 98%`)[[1]] + table(results_suppressed2$`significance 98%`)[[2]]), 3),
  round(table(results_suppressed2$`significance 98.5%`)[[1]] / (table(results_suppressed2$`significance 98.5%`)[[1]] + table(results_suppressed2$`significance 98.5%`)[[2]]), 3),
  round(table(results_suppressed2$`significance 99%`)[[1]] / (table(results_suppressed2$`significance 99%`)[[1]] + table(results_suppressed2$`significance 99%`)[[2]]), 3)
), level = c('Significance 95%', 'Significance 97.5%', 'Significance 98%', 'Significance 98.5%', 'Significance 99%'), Accuracy = rep(c('Specificity'), 5)
)

graph_results_suppressed2<-ggplot(results_suppressed2, aes(x=`P-value`)) +
  geom_histogram(color="black", fill="red", position="dodge")+
  geom_vline( aes(xintercept=0.05),
              linetype="dashed", size = 2)  +
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
# graph_results_suppressed2
jpeg('other_figures/compare_value_with_previous.jpeg', units = "in", width = 8, height = 6, res = 300)
graph_results_suppressed2
dev.off()

# jpeg('other_figures/p_value_distributions.jpeg', units = "in", width = 15, height = 13, res = 300)
# ggpubr::ggarrange(graph_results_suppressed2, graph_results_suppressed1, graph_compare_first_peak_value,
#           labels = c("A", "B", "C"),
#           ncol = 2, nrow = 2,
#           font.label = list(size = 20, color = "black")
#           )
# dev.off()
jpeg('other_figures/report_figures_3.jpeg', units = "in", width = 15, height = 13, res = 300)
ggpubr::ggarrange(graph_results_suppressed2, graph_results_suppressed1, graph_compare_first_peak_value,
                  labels = c("A", "B", "C"),
                  ncol = 2, nrow = 2,
                  font.label = list(size = 20, color = "black")
)
dev.off()

accuracy_dataset <- data.frame(rbind(sensitivity_value, specificity_value)) %>%
  arrange(level) %>%
  ggplot(aes(x = level, y = Accuracy.values, fill = Accuracy)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c('#999999','#E69F00'))  +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    axis.text.x = element_text(angle = 60, hjust = 1)#,
    # legend.position = "none"
  )
# accuracy_dataset
jpeg('other_figures/accuracy_barplot.jpeg', units = "in", width = 8, height = 6, res = 300)
accuracy_dataset
dev.off()

##########Plotting the Z-scores###########
z_score_data <- compare_first_peak_value %>%
  dplyr::select(id, value, z_stat) %>%
  mutate(`Group` = 'to unsuppressed') %>%
  bind_rows(results_suppressed2 %>% dplyr::select(id, value, z_stat) %>% mutate(`Group` = 'fully suppressed')) 
accuracy_dataset <- c()
threshold <- seq(0,3, 0.2)
for (i in 1:length(threshold)) {
    dat_set1 <- z_score_data %>%
      filter(Group == 'to unsuppressed') %>%
    mutate(x = ifelse(as.numeric(z_stat) > threshold[i], TRUE, FALSE))
    sensitivity_value <- round(table(dat_set1$x)[[2]] / (table(dat_set1$x)[[1]] + table(dat_set1$x)[[2]]), 3)
    
    dat_set2 <- z_score_data %>%
      filter(Group == 'fully suppressed') %>%
      mutate(x = ifelse(as.numeric(z_stat) > threshold[i], TRUE, FALSE))
    specificity_value <- round(table(dat_set2$x)[[1]] / (table(dat_set2$x)[[1]] + table(dat_set2$x)[[2]]), 3)
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
# accuracy_graph
# jpeg('other_figures/accuracy_linegraph.jpeg', units = "in", width = 8, height = 6, res = 300)
# accuracy_graph
# dev.off()
jpeg('other_figures/report_figure4.jpeg', units = "in", width = 8, height = 6, res = 300)
accuracy_graph
dev.off()