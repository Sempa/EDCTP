library(readxl)
library(tidyverse)
library(gtsummary)
library(dplyr)
library(ggplot2)
library(minpack.lm)
library(nlme)
library(splines)
library(cowplot)

# source('estimate_measurement_noise_dist.R')
#############################################
#' Estimating Sedia LAg Measurement noise
#############################################
sedia_generic <- read_csv("data/20180410-EP-LAgSedia-Generic.csv")

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
  rename(sedia_ODn = result...72) %>%
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
    model_data[counter, ] <- ODn_regression_function(
      dat = subset(dat, dat$subject_label_blinded == unique(dat$subject_label_blinded)[i])
    )
    # print(i)
  }
  # summary(model_data$slope_link_identity, na.rm = T)

  jpeg(paste(folder_name, "/", "Figure_link_identity", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)

  print(
    model_data %>%
      ggplot(aes(x = slope_link_identity)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, bins = 30, position = "identity") +
      scale_fill_manual(values = c("#69b3a2", "#404080")) +
      labs(fill = "") +
      geom_density(
        lwd = 1.2,
        linetype = 1,
        colour = 2
      ) +
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
      ggtitle(paste("histogram mean slope link=identity_", threshold))
  )
  dev.off()
  # summary(model_data$slope_link_log, na.rm = T)

  jpeg(paste(folder_name, "/", "Figure_link_log", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)

  print(
    model_data %>%
      ggplot(aes(x = slope_link_log)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, bins = 30, position = "identity") +
      scale_fill_manual(values = c("#69b3a2", "#404080")) +
      labs(fill = "") +
      geom_density(
        lwd = 1.2,
        linetype = 1,
        colour = 2
      ) +
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
      ggtitle(paste("histogram mean slope link=log_", threshold))
  )
  dev.off()

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
results <- c()
for (i in c(100, 400, 1000)) {
  x <- slopes_for_suppressed(ODn_vl_data = sedia_generic, threshold = i)
  results <- rbind(results, x)
}
write.csv(results, "output_table/results_suppressed.csv") # results
