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
  # filter(Group == 'early suppressions') %>%
  select(subject_label_blinded, days_since_eddi, `days since tx start`, Sex = BiologicalSex, 
         Age = `Age at Draw`, test_date, sedia_ODn, viral_load, visits, Group) # , baseline_visit 

full_dataset <- bind_rows(cephia_pts %>%
  mutate(cohort = 'cephia'),
  read.csv('data/full_africos_data_with_ODn.csv') %>% #africos_pts <- 
  mutate(days_since_eddi = NA,
         `days since tx start` = NA,
         Group = ifelse(exclude == 1, 'early suppressions', ifelse(is.na(exclude), 'suppression failures', NA)),
         Group = ifelse(is.na(Group), 'suppression failures', Group)) %>%
  select(subject_label_blinded = subjid, days_since_eddi, `days since tx start`, Sex, Age, 
         test_date = Date_Collected, sedia_ODn = ODn, viral_load = vl, 
         visits = Visit, Group, exclude) %>%
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
         subject_label_blinded = as.double(cur_group_id()),
         Age = as.numeric(Age))%>%
  select(subject_label_blinded, days_since_eddi, `days since tx start`, Sex, Age, 
         test_date, sedia_ODn, viral_load, visits, Group, cohort)
)


