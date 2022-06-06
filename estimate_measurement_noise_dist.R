library(readxl)
library(tidyverse)
library(gtsummary)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(inctools)
library(plotrix) ## helps to estimate std. error so as to compare the means between kit calibrators
library(Exact)
library(feather)
library(doRNG)
library(janitor)
library(grDevices)
library(Cairo)
library(ggpmisc)

sedia_generic <- read_csv("data/20180410-EP-LAgSedia-Generic.csv") # %>%
  # filter(visit_hivstatus == "P") %>%
  # rename(sedia_ODn = result...72) %>%
  # filter(visit_id != 21773 & visit_id != 21783 & visit_id != 21785)

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
  mutate(`sigma over mean ODn` = `sigma Sedia ODn`/`mean Sedia ODn`) 

sd_sedia_ODn_vs_ODn_bc <- sedia_distribution_blinded_controls %>%
  ggplot(aes(x = `mean Sedia ODn`, y = `sigma Sedia ODn`)) +
  geom_point() +
  expand_limits(x = 0, y = 0) +
  geom_smooth(method = lm, size = 1.5, se = FALSE)

sd_sedia_ODn_vs_ODn_bc

noise <- summary(glm(sd_Sedia_ODn ~ mean_Sedia_ODn, data = sedia_distribution_blinded_controls %>%
  mutate(sd_Sedia_ODn = `sigma Sedia ODn`, mean_Sedia_ODn = `mean Sedia ODn`)
  )
  )

sigma_ODn_vs_ODn_bc <- sedia_distribution_blinded_controls %>%
  ggplot(aes(x = `mean Sedia ODn`, y = `sigma over mean ODn`)) +
  geom_point() + 
  expand_limits(x = 0, y = 0) +
  geom_smooth(method = lm, size = 1.5, se = FALSE)

sigma_ODn_vs_ODn_bc
summary(glm(sigma_ODn ~ mean_Sedia_ODn, data = sedia_distribution_blinded_controls %>%
              mutate(sigma_ODn = `sigma over mean ODn`, mean_Sedia_ODn = `mean Sedia ODn`)
)
)


close_per_patient_visits <- sedia_generic %>%
  filter(visit_hivstatus == "P") %>%
  rename(sedia_ODn = result...72) %>%
  filter(visit_id != 21773 & visit_id != 21783 & visit_id != 21785) %>%
  select(subject_label_blinded, test_date, sedia_ODn, viral_load, art_initiation_date, aids_diagnosis_date, 
         art_interruption_date, art_resumption_date, treatment_naive, 
         on_treatment, first_treatment) %>%
  filter(!is.na(sedia_ODn)) %>%
  distinct(subject_label_blinded, test_date, .keep_all = T) %>%
  arrange(subject_label_blinded, test_date) %>%
  group_by(subject_label_blinded) %>%
  mutate(inter_test_interval = c(0, diff(test_date)),
         previous_date = test_date- inter_test_interval
         ) %>%
  mutate(inter_test_interval = c(diff(test_date),0),
         next_date = test_date+ inter_test_interval
  ) %>%
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
write.csv(close_per_patient_visits, 'data/close_per_patient_visits.csv')
close_per_patient_visits_new <- read_csv('data/close_per_patient_visits_new.csv') %>%
  mutate(id = paste(subject_label_blinded, code)) %>%
  group_by(id) %>%
  summarise(
    `mean Sedia ODn` = round(mean(sedia_ODn), 3),
    `sigma Sedia ODn` = round(sd(sedia_ODn), 3)
  ) %>%
  mutate(`sigma over mean ODn` = `sigma Sedia ODn`/`mean Sedia ODn`) 

sd_sedia_ODn_vs_ODn<- close_per_patient_visits_new %>%
  ggplot(aes(x = `mean Sedia ODn`, y = `sigma Sedia ODn`)) +
  geom_point() +
  expand_limits(x = 0, y = 0) +
  geom_smooth(method = lm, size = 1.5, se = FALSE)

sd_sedia_ODn_vs_ODn
summary(glm(sigma_Sedia_ODn ~ mean_Sedia_ODn, data = close_per_patient_visits_new %>%
              mutate(sigma_Sedia_ODn = `sigma Sedia ODn`, mean_Sedia_ODn = `mean Sedia ODn`)
)
)

sigma_ODn_vs_ODn <- close_per_patient_visits_new %>%
  ggplot(aes(x = `mean Sedia ODn`, y = `sigma over mean ODn`)) +
  geom_point() +
  expand_limits(x = 0, y = 0) +
  geom_smooth(method = lm, size = 1.5, se = FALSE)

sigma_ODn_vs_ODn
summary(glm(sigma_over_mean_ODn ~ mean_Sedia_ODn, data = close_per_patient_visits_new %>%
              mutate(sigma_over_mean_ODn = `sigma over mean ODn`, mean_Sedia_ODn = `mean Sedia ODn`)
            )
        )
