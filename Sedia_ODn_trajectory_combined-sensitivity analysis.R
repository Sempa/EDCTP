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

sedia_generic <- read_csv("data/20180410-EP-LAgSedia-Generic.csv") %>%
  mutate(sedia_ODn = `result...15`)

slope_data1 <- readRDS('data/model_data.rds') %>%
  mutate(subject_label_blinded = id) %>%
  left_join(sedia_generic %>%
              select(
                subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load, 
                sex, age_in_years, days_since_first_art, subtype
              ) %>%
              mutate(s_type = as.factor(ifelse(subtype == 'B', 'B', ifelse(subtype == 'C', 'C', 'Other')))),
            by = c('subject_label_blinded')) %>%
  distinct_at(vars(subject_label_blinded), .keep_all = T)
summary(lm(slope_link_identity ~ sex + age_in_years + s_type, data = slope_data1))
glm(slope_link_identity ~ sex + age_in_years + s_type, data = slope_data1) %>%tbl_regression()
anova(lm(slope_link_identity ~ sex + age_in_years + s_type, data = slope_data1))

slope_data2 <- readRDS('data/model_data_from_suppressed.rds') %>%
  mutate(subject_label_blinded = as.numeric(substring(id, 1,8))) %>%
  distinct_at(vars(subject_label_blinded), .keep_all = T) %>%
  left_join(sedia_generic %>%
              select(
                subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load, 
                sex, age_in_years, days_since_first_art, subtype
              ) %>%
              mutate(s_type = as.factor(ifelse(subtype == 'B', 'B', ifelse(subtype == 'C', 'C', 'Other')))),
            by = c('subject_label_blinded')) %>%
  distinct_at(vars(subject_label_blinded), .keep_all = T) %>%
  mutate(slope_link_identity = as.numeric(slope_link_identity))
summary(lm(slope_link_identity ~ sex + age_in_years + s_type, data = slope_data2))
glm(slope_link_identity ~ sex + age_in_years + s_type, 
    data = slope_data2, 
    family = gaussian(link = "identity")) %>%
  tbl_regression()
anova((lm(slope_link_identity ~ sex + age_in_years + s_type, data = slope_data2)))

# differences
sedia_data_for_diff <- readRDS('data/data_for_differences.rds') %>%
  left_join(sedia_generic %>%
              select(
                subject_label_blinded, test_date,
                sex, age_in_years, subtype
              ) %>%
              mutate(s_type = as.factor(ifelse(subtype == 'B', 'B', ifelse(subtype == 'C', 'C', 'Other')))),
            by = c('subject_label_blinded', 'test_date'))
summary(nlme::lme(sedia_diff ~ days_since_eddi + sex + age_in_years + s_type, 
                  random = ~ 1|subject_label_blinded,
                  data = sedia_data_for_diff %>% filter(suprressed_throughout_followup_1000==1)))
anova((nlme::lme(sedia_diff ~ days_since_eddi + sex + age_in_years + s_type, 
                 random = ~ 1|subject_label_blinded,
                 data = sedia_data_for_diff %>% filter(suprressed_throughout_followup_1000==1))))
## to preak sedia diff comparisons
summary(lm(sedia_diff ~ days_since_eddi + sex + age_in_years + s_type, 
        data = sedia_data_for_diff %>% filter(to_peak==1)))
anova(lm(sedia_diff ~ days_since_eddi + sex + age_in_years + s_type, 
         data = sedia_data_for_diff %>% filter(to_peak==1)))

# Analysis for those with visits straddling around ART start
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
  select(subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load, 
         visits, baseline_visit, Subtype, age = `Age at Draw`, BiologicalSex)
summary(lme4::lmer(sedia_ODn ~ days_since_eddi + age + (1|subject_label_blinded),
                    data = final_dataset %>%mutate(sex = as.factor(BiologicalSex ))))
lme4::lmer(sedia_ODn ~ days_since_eddi + age + (1|subject_label_blinded),
            data = final_dataset %>%mutate(sex = as.factor(BiologicalSex ))) %>%
  tbl_regression()
car::Anova(lme4::lmer(sedia_ODn ~ days_since_eddi + age + (1|subject_label_blinded),
                      data = final_dataset %>%mutate(sex = as.factor(BiologicalSex ))))