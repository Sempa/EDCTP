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
  # arrange(QC_9)
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
dt03 <- bind_cols(dt01, calibratorData_sample) %>%
  mutate(diff = OD-as.numeric(value))
  # full_join(dt02, by = 'id') %>%
  # arrange(`Study Number`, `Date Collected`)
dt04 <- africos_calibratorData %>%
  filter(calibrator_infor != 'sample') %>%
  group_by(Plate, calibrator_infor) %>%
  summarise(median_OD = median(as.numeric(value))) %>%
  ungroup() %>%
  filter(calibrator_infor == 'CAL')
dt05 <- dt03 %>%
  left_join(dt04, by = 'Plate')