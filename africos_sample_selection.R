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

pt_dataset <- read_excel("data/rv329_sacema_12sep22_Plasma Inventory[3].xlsx", sheet = "Request")[,-c(19:21)] %>%
  # mutate(id_visit = paste(`SUBJECT ID (CHAR)`, `STUDY VISIT (NUM)`, sep = '-')) %>%
  filter(!(`SUBJECT ID (CHAR)` == 'A01-0008' & `STUDY VISIT (NUM)` == 18 & `VL Copies/mL` == 50 & vl == 25)) %>%
  filter(!(`SUBJECT ID (CHAR)` == 'A01-0017' & `STUDY VISIT (NUM)` == 17 & `VL Copies/mL` == 1740 & vl == 1740)) %>%
  filter(!(`SUBJECT ID (CHAR)` == 'A01-0152' & `STUDY VISIT (NUM)` == 15 & `VL Copies/mL` == 50 & vl == 25)) %>%
  filter(!(`SUBJECT ID (CHAR)` == 'A01-0204' & `STUDY VISIT (NUM)` == 14 & is.na(`VL Copies/mL`) & vl == 1)) %>%
  filter(!(`SUBJECT ID (CHAR)` == 'A01-0224' & `STUDY VISIT (NUM)` == 15 & is.na(`VL Copies/mL`) & vl == 1)) %>%
  filter(!(`SUBJECT ID (CHAR)` == 'A01-0225' & `STUDY VISIT (NUM)` == 14 & `VL Copies/mL` == 88 & vl == 88)) %>%
  filter(!(`SUBJECT ID (CHAR)` == 'A01-0226' & `STUDY VISIT (NUM)` == 13 & `VL Copies/mL` == 102 & vl == 102)) %>%
  mutate(`STUDY VISIT (NUM)` == ifelse((`STUDY VISIT (NUM)` == 'A01-0020' & `VL Copies/mL` == 2116 & vl == 2116), 3, `STUDY VISIT (NUM)`))
write_rds(pt_dataset, 'data/africos_cohort.rds')
pt_samples <- read_excel("data/rv329_sacema_12sep22_Plasma Inventory[3].xlsx", sheet = "Sheet1")[-c(1:3),] %>%
  arrange(`Study Number`) %>%
  mutate(Uganda = grepl('A', `Study Number`),
         samples_exist = 1
         ) %>% #select samples from Uganda only for logistical purposes
  filter(Uganda == T) %>%
  dplyr::select(`Study Number`, Plasma_acd_volume, plasma_edta_volume, Plasma_acd_vailcount, plasma_edta_vailcount, samples_exist)

final_dataset <- pt_dataset %>%
  mutate(Uganda = grepl('A', `SUBJECT ID (CHAR)`))%>% #select samples from Uganda only for logistical purposes
  filter(Uganda == T) %>%
  group_by(`SUBJECT ID (CHAR)`) %>%
  mutate(art_start = min(as.Date(`NUMERIC DATE FOR art_sdt`, '%d/%m/%y'), na.rm = T),
         first_visit_date = min(as.Date(`Visit Date`, '%d/%m/%y')),
         flag_visit = ifelse(!is.na(as.Date(`NUMERIC DATE FOR art_sdt`, '%d/%m/%y')), `STUDY VISIT (NUM)`, NA) #visit_first_recorded_art = 
         ) %>%
  mutate(min_flag_visit = min(flag_visit, na.rm = T),
         visit_after_ARTstart = ifelse(`STUDY VISIT (NUM)` > min_flag_visit, 1, 0)) %>%
  mutate(duration_on_ART = max(`Duration of started ART (years)`, na.rm = T),
         visits_before_art = ifelse(art_start>`Visit Date`, 1,0),
         ) %>%
  ungroup() %>%
  mutate(stradlling_art_start = ifelse(`NUMERIC DATE FOR art_sdt` >= first_visit_date, T, F), #166 with ART start after first VL visit
         `missing vl` = ifelse(!is.na(vl), F, T), # 236 missing viral loads But 38 visits have plasma vl available
         `morethan 2yrs` = ifelse(duration_on_ART >=2, 1, 0) ## All have been followed up for at least 2 years
         ) %>% 
  # filter(!is.na(vl)) %>%
  mutate(detectable_vl = ifelse((visit_after_ARTstart == 1 & `VL Copies/mL` >1000), 1, 0)) %>%
  group_by(`SUBJECT ID (CHAR)`) %>%
  mutate(ever_unsuppressed_after_artstart = max(detectable_vl, na.rm = T),
         `Study Number` = paste(`SUBJECT ID (CHAR)`, `STUDY VISIT (NUM)`, sep = '-')
         ) %>% #28 people with detectable VL after ART start
  right_join(pt_samples, by = 'Study Number') %>%
  filter(samples_exist==1)

# m=final_dataset%>%dplyr::select(`SUBJECT ID (CHAR)`, `PlACD Vol`, `PlACD  Vials` , vl, `missing vl`) %>% mutate(`missing plasma` = ifelse((`PlACD Vol`!='N/A' | `PlACD  Vials` != 'N/A') & is.na(vl), 1,0))
# table(m$`missing plasma`)
# x=final_dataset%>%
#    filter(ever_unsuppressed_after_artstart == 1) #%>%
# length(unique((final_dataset%>% mutate(country = grepl('A', `SUBJECT ID (CHAR)`)) %>% filter(country ==T & ever_unsuppressed_after_artstart==1))$`SUBJECT ID (CHAR)`))
x_data <- pt_dataset %>%
  filter(`SUBJECT ID (CHAR)` %in% unique((final_dataset%>% 
                                            mutate(country = grepl('A', `SUBJECT ID (CHAR)`)) %>% 
                                            filter(ever_unsuppressed_after_artstart==1))$`SUBJECT ID (CHAR)`))
write.csv(x_data, 'data/ever_unsuppressed.csv')
history_deteable_vl <- read.csv('data/ever_unsuppressed_edited.csv')
samples_to_be_test <- history_deteable_vl %>%
  filter(!is.na(visit_to.selct_from)) %>%
  dplyr::select(`PlACD Vol` = PlACD.Vol, `PlACD Vials` = PlACD..Vials, `PLEDTA.Vol` = PLEDTA.Vol, 
                `PLEDTA Vials` = PLEDTA..Vials, `SUBJECT ID` = SUBJECT.ID..CHAR.,
                `Visit Date` = Visit.Date, `STUDY VISIT NUM` = STUDY.VISIT..NUM., 
                `NUMERIC DATE FOR art_sdt` = NUMERIC.DATE.FOR.art_sdt, 
                `Duration of started ART years` = Duration.of.started.ART..years., 
                `HIV status at visit` = HIV.status.at.visit,
                `Viral Load Sample Draw Period` = Viral.Load.Sample.Draw.Period, 
                `Viral Load Result Prefix` = Viral.Load.Result.Prefix, `VL Copies mL` = VL.Copies.mL, vl)
write.csv(samples_to_be_test, 'data/List of samples to test.csv')
