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
  select(id, to_peak, set_to_peak, to_trough, set_to_trough)
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
    to_peak, set_to_peak, to_trough, set_to_trough
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

table_1 <- table(sedia_slope_data$slope_cat, sedia_slope_data$type)

A <- table_1[1,2]; B <- table_1[1,1]
C <- table_1[2,2]; D <- table_1[2,1]
             
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
  
  table_1 <- table(sedia_slope_data$slope_cat, sedia_slope_data$type)
  
  A <- table_1[1,2]; B <- table_1[1,1]
  C <- table_1[2,2]; D <- table_1[2,1]
  
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

# model1 <- glm(slope_cat ~ type, 
#          data = sedia_slope_data %>%
#            mutate(type = as.factor(flag_viralload_detectablity2)), 
#          family = binomial())
# threshold=0.5
# predicted_values<-ifelse(predict(model1, type="response")>threshold,1,0)
# actual_values<-model1$y
# conf_matrix<-table(predicted_values,actual_values)
# conf_matrix
# caret::sensitivity(conf_matrix)

caret::confusionMatrix(table(sedia_slope_data$slope_cat, sedia_slope_data$type))

  mutate(xi= test_date - previous_date,
         xii= next_date - test_date,
         min_date = min(test_date)) %>%
  mutate(flag = ifelse(xi<=10, 1, ifelse(xii<=10 , 1, 0))) %>% #& min_date!=test_date & min_date!=test_date
  mutate(flag2 = ifelse(xi>10 & xii<=10,1,0)) %>%
  # filter(flag == 1) %>%
  mutate(visits = 1:length(subject_label_blinded)) %>%
  mutate(n_visits = max(visits)) %>%
  filter(n_visits>1) %>%
  select(subject_label_blinded, test_date, xi, previous_date, xii, next_date, flag, flag2, sedia_ODn, viral_load, art_initiation_date, aids_diagnosis_date, 
         art_interruption_date, art_resumption_date, treatment_naive, 
         on_treatment, first_treatment) %>%
  ungroup() %>%
  mutate(y1= ifelse(xi<=10 & xi>0, 1,0),
         y2= ifelse(xii<=10 & xii>0, 1,0)) %>%
  select(subject_label_blinded, test_date, xi, previous_date, xii, next_date, flag, flag2, sedia_ODn, viral_load, y1, y2) %>%
  filter(y1==1 | y2==1) %>%
  group_by(subject_label_blinded) %>%
  mutate(inter_test_interval = c(0, diff(test_date)),
         previous_date = test_date- inter_test_interval
  ) %>%
  mutate(inter_test_interval = c(diff(test_date),0),
         next_date = test_date+ inter_test_interval
  ) %>%
  mutate(xi= test_date - previous_date,
         xii= next_date - test_date,
         min_date = min(test_date))