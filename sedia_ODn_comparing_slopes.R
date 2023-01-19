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
sedia_generic <- read_csv("data/20180410-EP-LAgSedia-Generic.csv")# %>%
  # mutate(sedia_ODn = `result...15`)
data_intermitent_suppression_selected_visits <- read.csv("output_table/intermitent_suppression_selected.csv") %>%
  filter(!is.na(to_peak) |!is.na(to_trough)) %>%
  mutate(id = paste(subject_label_blinded, days_since_eddi, sep = '_'),
         id_2 = id) %>%
  select(id, to_peak, set_to_peak, to_peak_first_visit, to_trough, set_to_trough, to_trough_first_visit)
data_sorting <- sedia_generic %>%
  filter(visit_hivstatus == "P") %>%
  rename(sedia_ODn = result...72) %>%
  filter(visit_id != 21773 & visit_id != 21783 & visit_id != 21785) %>%
  select(
    subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load, art_initiation_date, aids_diagnosis_date,
    art_interruption_date, art_resumption_date, treatment_naive,
    on_treatment, first_treatment
  ) %>%
  filter(!is.na(sedia_ODn)) %>%
  mutate(id = paste(subject_label_blinded, days_since_eddi, sep = "_")) %>%
  full_join(data_intermitent_suppression_selected_visits, by = "id") %>%
  arrange(subject_label_blinded, days_since_eddi) %>%
  group_by(subject_label_blinded) %>%
  mutate(
    flagvl_100 = ifelse(viral_load <= 100, 0, 1),
    flagvl_400 = ifelse(viral_load <= 400, 0, 1),
    flagvl_1000 = ifelse(viral_load <= 1000, 0, 1)
  ) %>%
  mutate(
    mean_flag_100 = mean(flagvl_100),
    mean_flag_400 = mean(flagvl_400),
    mean_flag_1000 = mean(flagvl_1000),
  ) %>%
  mutate(
    suprressed_throughout_followup_100 = ifelse(mean_flag_100 == 0, 1, 0),
    suprressed_throughout_followup_400 = ifelse(mean_flag_400 == 0, 1, 0),
    suprressed_throughout_followup_1000 = ifelse(mean_flag_1000 == 0, 1, 0)
  ) %>%
  mutate(
    eddi_diff = c(diff(days_since_eddi), 0),
    days_since_eddi_2 = days_since_eddi + eddi_diff
  ) %>%
  mutate(
    sedia_ODn_diff = c(diff(sedia_ODn), 0),
    sedia_ODn_2 = sedia_ODn + sedia_ODn_diff
  ) %>%
  mutate(
    viraload_diff = c(diff(viral_load), 0),
    viral_load_2 = viral_load + viraload_diff
  ) %>%
  # mutate(eddi = )
  select(
    subject_label_blinded, days_since_eddi, days_since_eddi_2, sedia_ODn, sedia_ODn_2, test_date, viral_load, viral_load_2,
    suprressed_throughout_followup_100, suprressed_throughout_followup_400, suprressed_throughout_followup_1000,
    to_peak, set_to_peak, to_peak_first_visit, to_trough, set_to_trough, to_trough_first_visit
  ) %>%
  filter(days_since_eddi != days_since_eddi_2) %>%
  arrange(subject_label_blinded, days_since_eddi) %>%
  mutate(flag_viralload_detectablity1 = ifelse(viral_load_2 <= 1000 | viral_load <= 1000, 1, 0)) %>%
  mutate(flag_viralload_detectablity2 = ifelse(viral_load <= 1000 & viral_load_2 <= 1000, "suppressed", 
                                               ifelse(viral_load <= 1000 & viral_load_2 > 1000, "to_fro_peak", 
                                                      ifelse(viral_load > 1000 & viral_load_2 <= 1000, "to_fro_peak", NA)
                                                      )
                                               )
         ) %>%
  mutate(sedia_slope = ifelse(flag_viralload_detectablity2 != 'NA', (sedia_ODn_2 - sedia_ODn)/(days_since_eddi_2 - days_since_eddi)), NA)

accuracy_function_all <- function(dat, threshold){
sedia_slope_data <- dat %>%
  filter(!is.na(sedia_slope)) %>%
  mutate(slope_cat = ifelse(sedia_slope<= threshold, 1, 0),
         type = as.factor(flag_viralload_detectablity2))

table_1 <- table(sedia_slope_data$type, sedia_slope_data$slope_cat)

A <- table_1[2,2]; B <- table_1[2,1]
C <- table_1[1,2]; D <- table_1[1,1]
             
sensitivity <- A/(A+C)

specificity <- D/(B+D)

Prevalence <- (A+C)/(A+B+C+D)

PPV <- (sensitivity * Prevalence)/((sensitivity*Prevalence) + ((1-specificity)*(1-Prevalence)))

NPV <- (specificity * (1-Prevalence))/(((1-sensitivity)*Prevalence) + ((specificity)*(1-Prevalence)))

return( cbind(threshold = threshold,
              sensitivity = sensitivity,
              specificity = specificity,
              # Prevalence = Prevalence,
              PPV = PPV,
              NPV = NPV)
        )
}

accuracy_data_all <- c()
counter <- 0
threshold <- c(0.0001, 0.001, 0.002)
for (i in 1:length(threshold)) {
  results <- accuracy_function_all(dat = data_sorting,
                               threshold = threshold[i])
  accuracy_data_all <- rbind(accuracy_data_all, results)
}
accuracy_data_all

sedia_slope_data_2 <- data_sorting %>%
  filter(suprressed_throughout_followup_1000 == 1 | to_peak == 1 | to_trough == 1
         ) %>%
  filter(!is.na(sedia_slope))

accuracy_function <- function(dat, threshold){
  sedia_slope_data <- dat %>%
    filter(suprressed_throughout_followup_1000 == 1 | to_peak == 1 | to_trough == 1
    ) %>%
    filter(!is.na(sedia_slope)) %>%
    mutate(slope_cat = ifelse(sedia_slope<= threshold, 1, 0),
           type = as.factor(flag_viralload_detectablity2))
  
  table_1 <- table(sedia_slope_data$type, sedia_slope_data$slope_cat)
  
  A <- table_1[2,2]; B <- table_1[2,1]
  C <- table_1[1,2]; D <- table_1[1,1]
  
  sensitivity <- A/(A+C)
  
  specificity <- D/(B+D)
  
  Prevalence <- (A+C)/(A+B+C+D)
  
  PPV <- (sensitivity * Prevalence)/((sensitivity*Prevalence) + ((1-specificity)*(1-Prevalence)))
  
  NPV <- (specificity * (1-Prevalence))/(((1-sensitivity)*Prevalence) + ((specificity)*(1-Prevalence)))
  
  return( cbind(threshold = threshold,
                sensitivity = sensitivity,
                specificity = specificity,
                # Prevalence = Prevalence,
                PPV = PPV,
                NPV = NPV)
  )
}

accuracy_data <- c()
counter <- 0
threshold <- c(0.0001, 0.001, 0.002)
for (i in 1:length(threshold)) {
  results <- accuracy_function(dat = data_sorting,
                                   threshold = threshold[i])
  accuracy_data <- rbind(accuracy_data, results)
}
accuracy_data
write.csv(accuracy_data, 'output_table/accuracy_data.csv')
######################################################
#' The differences between Sedia LAg ODn visits
######################################################
sedia_diffs_data_sorting <- data_sorting %>%
  filter(suprressed_throughout_followup_1000 == 1 | to_peak == 1 | to_trough == 1
  ) %>%
  filter(!is.na(sedia_slope)) %>%
  arrange() %>%
  mutate(sedia_diff = sedia_ODn_2 - sedia_ODn)

t.test((sedia_diffs_data_sorting %>%
          filter(suprressed_throughout_followup_100 == 1))$sedia_diff)
t.test((sedia_diffs_data_sorting %>%
          filter(suprressed_throughout_followup_400 == 1))$sedia_diff)
t.test((sedia_diffs_data_sorting %>%
          filter(suprressed_throughout_followup_1000 == 1))$sedia_diff)
t.test((sedia_diffs_data_sorting %>%
         filter(to_peak_first_visit == 1))$sedia_diff)
t.test((sedia_diffs_data_sorting %>%
         filter(to_trough_first_visit == 1))$sedia_diff)

dat <- sedia_diffs_data_sorting %>%
  filter(suprressed_throughout_followup_100 == 1)
hist_plot1 <- hist(dat$sedia_diff, breaks = 30, freq = FALSE, ylab = '', xlab = 'Sedia LAg differences', main = '', col = 'black', cex.lab = 1.5, cex.axis = 1.5) # , main = paste("Mean slope link=identity_", threshold)
area_hist1 <- sum(hist_plot1$counts* abs(hist_plot1$mids[1]-hist_plot1$mids[2]))
n_values1 <- sum(hist_plot1$counts / area_hist1)
set.seed(11)
x <- rnorm(n_values1, mean = mean(dat$sedia_diff, na.rm = T), sd = sd(dat$sedia_diff))
dx <- density(x)
jpeg('other_figures/Figure_less_100.jpeg', units = "in", width = 8, height = 6, res = 300)
hist(dat$sedia_diff, breaks = 30, freq = FALSE, ylab = '', xlab = 'Sedia LAg differences', main = '', col = 'black', cex.lab = 2, cex.axis = 2) # , main = paste("Density plot link=identity_", threshold)
lines(dx, lwd = 3, col = "red")
dev.off()
  
dat <- sedia_diffs_data_sorting %>%
  filter(suprressed_throughout_followup_400 == 1)
hist_plot2 <- hist(dat$sedia_diff, breaks = 30, freq = FALSE, ylab = '', xlab = 'Sedia LAg differences', main = '', col = 'black', cex.lab = 1.5, cex.axis = 1.5) # , main = paste("Mean slope link=identity_", threshold)
area_hist2 <- sum(hist_plot2$counts* abs(hist_plot2$mids[1]-hist_plot2$mids[2]))
n_values2 <- sum(hist_plot2$counts / area_hist2)
set.seed(11)
x <- rnorm(n_values2, mean = mean(dat$sedia_diff, na.rm = T), sd = sd(dat$sedia_diff))
dx <- density(x)
jpeg('other_figures/Figure_less_400.jpeg', units = "in", width = 8, height = 6, res = 300)
hist(dat$sedia_diff, breaks = 30, freq = FALSE, ylab = '', xlab = 'Sedia LAg differences', main = '', col = 'black', cex.lab = 2, cex.axis = 2) # , main = paste("Density plot link=identity_", threshold)
lines(dx, lwd = 3, col = "red")
dev.off()
  
dat1 <- sedia_diffs_data_sorting %>%
  filter(suprressed_throughout_followup_1000 == 1)
hist_plot3 <- hist(dat1$sedia_diff, breaks = 30, freq = FALSE, ylab = '', xlab = 'Sedia LAg differences', main = '', col = 'black', cex.lab = 1.5, cex.axis = 1.5) # , main = paste("Mean slope link=identity_", threshold)
area_hist3 <- sum(hist_plot3$counts* abs(hist_plot3$mids[1]-hist_plot3$mids[2]))
n_values3 <- sum(hist_plot3$counts / area_hist3)
set.seed(11)
x1 <- rnorm(n_values3, mean = mean(dat1$sedia_diff, na.rm = T), sd = sd(dat1$sedia_diff))
dx1 <- density(x1)
# jpeg('other_figures/Figure_less_1000.jpeg', units = "in", width = 8, height = 6, res = 300)
# hist(dat1$sedia_diff, breaks = 30, freq = FALSE, xlim = c(-1, 2), ylab = '', xlab = 'Sedia LAg differences', main = '', col = 'black', cex.lab = 2, cex.axis = 2) # , main = paste("Density plot link=identity_", threshold)
# lines(dx1, lwd = 3, col = "red")
# dev.off()

dat2 <- sedia_diffs_data_sorting %>%
  filter(to_peak_first_visit == 1)
hist_plot4 <- hist(dat2$sedia_diff, breaks = 30, freq = FALSE, ylab = '', xlab = 'Sedia LAg differences', main = '', col = 'black', cex.lab = 1.5, cex.axis = 1.5) # , main = paste("Mean slope link=identity_", threshold)
area_hist4 <- sum(hist_plot4$counts* abs(hist_plot4$mids[1]-hist_plot4$mids[2]))
n_values4 <- sum(hist_plot4$counts / area_hist4)
set.seed(11)
x2 <- rnorm(n_values4, mean = mean(dat2$sedia_diff, na.rm = T), sd = sd(dat2$sedia_diff))
dx2 <- density(x2)
# jpeg('other_figures/Figure_to_peak.jpeg', units = "in", width = 8, height = 6, res = 300)
# hist(dat2$sedia_diff, breaks = 30, freq = FALSE, ylab = '', xlab = 'Sedia LAg differences', main = '', col = 'black', cex.lab = 2, cex.axis = 2) # , main = paste("Density plot link=identity_", threshold)
# lines(dx2, lwd = 3, col = "red")
# dev.off()
  
dat3 <- sedia_diffs_data_sorting %>%
  filter(to_trough_first_visit == 1)
hist_plot5 <- hist(dat3$sedia_diff, breaks = 30, freq = FALSE, ylab = '', xlab = 'Sedia LAg differences', main = '', col = 'black', cex.lab = 1.5, cex.axis = 1.5) # , main = paste("Mean slope link=identity_", threshold)
area_hist5 <- sum(hist_plot5$counts* abs(hist_plot5$mids[1]-hist_plot5$mids[2]))
n_values5 <- sum(hist_plot5$counts / area_hist5)
set.seed(11)
x3 <- rnorm(n_values5, mean = mean(dat3$sedia_diff, na.rm = T), sd = sd(dat3$sedia_diff))
dx3 <- density(x3)

jpeg('other_figures/Figure_combined_3.jpeg', units = "in", width = 15, height = 13, res = 300)

par(mfrow = c(2, 2))

hist(dat1$sedia_diff, breaks = 30, freq = FALSE, xlim = c(-1, 2), ylab = '', xlab = 'Sedia LAg differences', main = 'A', col = 'black', cex.lab = 2, cex.axis = 2, cex.main = 2) # , main = paste("Density plot link=identity_", threshold)
lines(dx1, lwd = 3, col = "red")


hist(dat2$sedia_diff, breaks = 30, freq = FALSE, ylab = '', xlab = 'Sedia LAg differences', main = 'B', col = 'black', cex.lab = 2, cex.axis = 2, cex.main = 2) # , main = paste("Density plot link=identity_", threshold)
lines(dx2, lwd = 3, col = "red")

hist(dat3$sedia_diff, breaks = 30, freq = FALSE, ylab = '', xlab = 'Sedia LAg differences', main = 'C', col = 'black', cex.lab = 2, cex.axis = 2, cex.main = 2) # , main = paste("Density plot link=identity_", threshold)
lines(dx3, lwd = 3, col = "red")

dev.off()
  