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

africos_data <- readxl::read_excel("data/Sedia® HIV-1 LAg-Avidity EIA ODs results 2.xlsx") %>%
  filter(!is.na(`Date Collected`)) %>%
  mutate(id = paste(`Date Collected`, `Study Number`, sep = '_'),
         `ODs 1` = `ODs`,
         `ODs 2` = `...13`,
         `ODs 3` = `...14`,
         QC_8 = QC...8,
         QC_9 = 1:length(`Date Collected`)) %>%
  dplyr::select(id, `Date Collected`, `Study Number`, QC_8, QC_9, Plates, `ODs 1`, `ODs 2`,`ODs 3`)# %>%
africos_data_1 <- readxl::read_excel("data/Sedia® HIV-1 LAg-Avidity EIA ODs results 2 (3).xlsx")[-1, ] %>%
  dplyr::select(ODs_new = ODs)
dt01 <- africos_data %>%
  pivot_longer(cols = c('ODs 1', 'ODs 2','ODs 3'),
               names_to = 'ODs',
               values_to = 'OD')
dt02 <- readxl::read_excel("data/AFRICOS_Sempa Data Pull_24Jul23.xlsx", 
                           sheet = "Sheet1") %>%
  mutate(id = paste(Date_Collected, subjid, sep = '_')) %>%
  group_by(subjid) %>%
  mutate(switch_date = as.Date(as.numeric(ifelse(`ARV Start Date` == '.', '', `ARV Start Date`)), 
                             origin = "1900-01-01"),
         art_start = as.Date(min(as.numeric(ifelse(`ARV Start Date` == '.', '', `ARV Start Date`)), na.rm = T), 
                             origin = "1900-01-01")) %>%
  ungroup() %>%
  dplyr::select(id, subjid, Visit, Age, Sex, `Viral Load`, Date_Collected, switch_date, art_start)
africos_calibratorData <- bind_rows(readxl::read_excel("data/Sedia® HIV-1 LAg-Avidity EIA ODs results 2 calibrator data.xlsx", 
                                             sheet = "Sheet1") %>% 
                                      mutate(sheet = 1),
                                    readxl::read_excel("data/Sedia® HIV-1 LAg-Avidity EIA ODs results 2 calibrator data.xlsx", 
                                               sheet = "Sheet2") %>% 
                                      mutate(sheet = 2),
                                    readxl::read_excel("data/Sedia® HIV-1 LAg-Avidity EIA ODs results 2 calibrator data.xlsx", 
                                               sheet = "Sheet3") %>% 
                                      mutate(sheet = 3),
                                    readxl::read_excel("data/Sedia® HIV-1 LAg-Avidity EIA ODs results 2 calibrator data.xlsx", 
                                               sheet = "Sheet4") %>% 
                                      mutate(sheet = 4),
                                    readxl::read_excel("data/Sedia® HIV-1 LAg-Avidity EIA ODs results 2 calibrator data.xlsx", 
                                               sheet = "Sheet5") %>% 
                                      mutate(sheet = 5)
) %>%
  filter(!(row == 'H' & column == 12))
calibratorData_sample <- africos_calibratorData %>%
  filter(calibrator_infor == 'sample')
dt03 <- bind_cols(dt01, africos_data_1) %>%
  group_by(id) %>%
  mutate(OD = median(OD)) %>%
  dplyr::select(id, `Date Collected`, `Study Number`, QC_8, QC_9, Plate = Plates, OD) %>%
  distinct(id, .keep_all = T)

dt04 <- africos_calibratorData %>%
  filter(calibrator_infor == 'CAL') %>%
  group_by(Plate) %>%
  summarise(calibrator_OD = median(as.numeric(value))) %>%
  ungroup()
dt05 <- dt03 %>%
  left_join(dt04, by = 'Plate') %>%
  mutate(visit_date = format(as.Date(`Date Collected`, format = "%Y-%m-%d"), "%d/%m/%Y"),
         ODn = OD/calibrator_OD,
         id = paste(visit_date, `Study Number`, sep = '_'))
dt06 <- read.csv('data/ever_unsuppressed_edited.csv') %>%
  mutate(visit_date = format(as.Date(Visit.Date, format = "%d/%m/%Y"), "%d/%m/%Y")) %>%
  mutate(id = paste(visit_date, SUBJECT.ID..CHAR., sep = '_'),
         visit_date_new = as.Date(Visit.Date, format = "%d/%m/%Y")) %>%
  right_join(dt05, by = 'id') %>%
  arrange(SUBJECT.ID..CHAR., visit_date_new)
# write.csv(dt06, 'data/ever_suppressed_with_ODns.csv')
dt06_1 <- read_delim("data/ever_suppressed_with_ODns.csv", 
                     delim = ";", escape_double = FALSE, trim_ws = TRUE)

dt07 <- read_excel("data/AFRICOS_Sempa Data Pull_24Jul23.xlsx") %>%
  dplyr::select(Date_Collected, subjid, Visit, Age, Sex, `Viral Load`, `ARV Class`, ARV, `ARV Start Date`) %>%
  mutate(visit_date = format(as.Date(Date_Collected, format = "%Y-%m-%d"), "%d/%m/%Y"),
         id = paste(visit_date, subjid, sep = '_')) #%>%
dt08 <- dt07 %>%
  full_join(dt06_1, by = 'id') %>%
  dplyr::select(Date_Collected, subjid, Visit, Age, Sex, `Viral Load`, 
                VL.Copies.mL, vl, flag, vftot, slow_VL_decline, vl_before_detect, 
                group_vl_before_detect, vl_after_detect, to_detect, group_to_detect, 
                to_undetect, group_to_undetect, visit_date.y, visit_date_new, ODn = ODn_,
                exclude, comment)
write.csv(dt08, 'data/full_africos_data_with_ODn.csv')
dt09 <- dt08 %>%
  filter(is.na(exclude)) %>% # removing IDs/records where plasma samples were unavailable for testing.
  filter(vl_before_detect == 1) %>%
  group_by(subjid) %>%
  mutate(n_visits = length(subjid),
         subject_label_blinded = cur_group_id()) %>%
  ungroup() %>%
  mutate(test_date = as.character(Date_Collected)) %>%
  rename(sedia_ODn = ODn,
         viral_load = `Viral Load`, visits = Visit,
         to_peak = to_detect, set_to_peak = group_to_detect, to_trough = to_undetect, 
         set_to_trough = group_to_undetect, all_visits_to_peak = vl_before_detect) %>%
  dplyr::select(subject_label_blinded, test_date, sedia_ODn,
                viral_load, visits,
                to_peak, set_to_peak, to_trough, 
                set_to_trough, all_visits_to_peak)
write.csv(dt09, 'data/africos_data_with_ODn.csv')
