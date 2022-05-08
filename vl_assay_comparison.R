library(readxl)
library(tidyverse)
library(gtsummary)
library(dplyr)
library(ggplot2)
library(minpack.lm)
library(nlme)
library(splines)
library(cowplot)

final_dataset <- read_csv("Sempa_final_pull_with_results.csv") %>%
  mutate(vl = ifelse(`Viral Load at Draw`== '<40', '40', ifelse(`Viral Load at Draw`=='NULL', "", `Viral Load at Draw`))) %>%
  mutate(logvl = log(as.numeric(vl),10),
         time_on_trt = -1 * `days since tx start`,
         hiv_antibody_value = `Sedia LAg Odn screen`,
         RaceEthnicity = `Race/Ethnicity`) %>%
  arrange(subject_label_blinded, time_on_trt)

dt01 <- final_dataset %>%
  group_by(subject_label_blinded) %>%
  mutate(days_since_ART_start = max(time_on_trt, na.rm = T)) %>%
  ungroup() %>%
  distinct(subject_label_blinded, .keep_all = T) %>%
  dplyr::select(Group, Subtype, BiologicalSex, RaceEthnicity, days_since_ART_start) %>%
  tbl_summary(missing = "no") %>% 
  modify_header(label = "Variable") %>%
  bold_labels()
dt01


x=final_dataset[final_dataset$Group == 'early suppressions',] %>% 
  select(Group, subject_label_blinded, `days since tx start`, 
         `Viral Load at Draw`, vl, logvl, `Designated Elite at Draw`) %>% 
  arrange(subject_label_blinded, `days since tx start`)

plot_grid(
ggplot(data = final_dataset[final_dataset$Group == 'early suppressions',],
       aes(x = time_on_trt, y = logvl)) +
  geom_line(aes(color = as.factor(subject_label_blinded), linetype = Group), size = 1.5) +#aes(color = subject_label_blinded, linetype = Group), group = Group
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin=unit(c(0,0,0,0), "null"),
    legend.position = "none"
  ) +
  ylab(expression(paste(Log[10], ' Viral load'))) + xlab(''),

ggplot(data = final_dataset[final_dataset$Group == 'early suppressions',],
       aes(x = time_on_trt, y = `Sedia LAg Odn screen`)) + #as.factor(subject_label_blinded)
  geom_line(aes(color = as.factor(subject_label_blinded), linetype = Group), size = 1.5) +#aes(color = subject_label_blinded, linetype = Group), group = Group
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin=unit(c(0,0,0,0), "null"),
    legend.position = "none"
  ) +
  ylab('Sedia LAg ODn') + xlab('Time from ART start (Days)'),
labels = 'AUTO',
ncol = 1,
label_colour = "red"
)  
# graphing the LAg results straddling ART start (the so-called `baseline` or last 
# LAg reading before ART and the first LAg value at least 60months post ART start)
baseline_first_reading <- final_dataset %>%
  filter(time_on_trt <= 0) %>%
  group_by(subject_label_blinded) %>%
  mutate(last_LAg_reading_b4_ART = ifelse(time_on_trt <= 0, max(time_on_trt), 0),
         ) %>%
  filter(last_LAg_reading_b4_ART == time_on_trt) %>%
  dplyr::select(subject_label_blinded, time_on_trt) #, first_LAg_reading_on_ART, flag_visit

first_LAg_reading <- final_dataset %>%
  filter(time_on_trt>180.5) %>%
  group_by(subject_label_blinded) %>%
  mutate(
         first_LAg_reading_on_ART = ifelse(time_on_trt > 0, min(time_on_trt), 0)
  ) %>%
  filter(first_LAg_reading_on_ART == time_on_trt) %>%
  dplyr::select(subject_label_blinded, time_on_trt)

visits_stradling_ART_start <- rbind(baseline_first_reading, first_LAg_reading) %>%
  arrange(subject_label_blinded, time_on_trt)

LAg_stradling_ART_start <- final_dataset %>%
  dplyr::select(subject_label_blinded, time_on_trt, hiv_antibody_value, `Sedia LAg Odn screen`, logvl) %>%
  inner_join(visits_stradling_ART_start, by = c('subject_label_blinded', 'time_on_trt'))

# graphing LAg reading straddling ART start

ggplot(data = LAg_stradling_ART_start,
       aes(x = time_on_trt, y = `Sedia LAg Odn screen`)) + #as.factor(subject_label_blinded)
  geom_line(aes(color = as.factor(subject_label_blinded)), size = 1.5) +#aes(color = subject_label_blinded, linetype = Group), group = Group
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin=unit(c(0,0,0,0), "null"),
    legend.position = "none"
  ) +
  ylab('Sedia LAg ODn') + xlab('Time from ART start (Days)')

# studying the trajectory using CEPHIA samples that don't have pre 

assay_dataset <- read.csv("data/cephia_public_use_dataset_20210604.csv")

sedia_generic <- read_csv("data/20180410-EP-LAgSedia-Generic.csv") %>%
  filter(visit_hivstatus == "P") %>%
  rename(sedia_ODn = result...72) %>%
  filter(visit_id != 21773 & visit_id != 21783 & visit_id != 21785)
data_generated <- sedia_generic %>%
  select(subject_label_blinded, test_date, art_initiation_date, aids_diagnosis_date, 
         art_interruption_date, art_resumption_date, treatment_naive, 
         on_treatment, first_treatment, sedia_ODn, viral_load) %>%
  filter(!is.na(art_initiation_date)) %>%
  arrange(subject_label_blinded, test_date) %>%
  filter(!is.na(sedia_ODn)) %>%
  distinct(subject_label_blinded, test_date, .keep_all = T) %>%
  group_by(subject_label_blinded) %>%
  mutate(visits = 1:length(subject_label_blinded)) %>%
  mutate(n_visits = max(visits)) %>%
  filter(n_visits>1) %>%
  ungroup() %>%
  mutate(time_on_trt = as.numeric(test_date - art_initiation_date),
         `Sedia LAg Odn screen` = sedia_ODn) %>%
  dplyr::select(subject_label_blinded, test_date, art_initiation_date, time_on_trt, `Sedia LAg Odn screen`, viral_load, n_visits) %>%
  arrange(subject_label_blinded, test_date) 

set.seed(11)
  ggplot(data = data_generated %>%
           filter(subject_label_blinded %in% c(11018648, 11428415)), #%in% sample(unique(subject_label_blinded), size = 30, replace = F)
         aes(x = time_on_trt, y = `Sedia LAg Odn screen`)) + #as.factor(subject_label_blinded)
    geom_line(aes(color = as.factor(subject_label_blinded)), size = 1.5) +#aes(color = subject_label_blinded, linetype = Group), group = Group
    theme(
      text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.margin=unit(c(0,0,0,0), "null"),
      legend.position = "none"
    ) +
    ylab('Sedia LAg ODn') + xlab('Time from ART start (Days)')
  
# data_generated <- assay_dataset %>%
#   filter(assay == "LAg-Sedia",
#          assay_result_field == "final_result",
#          cephia_panel == "CEPHIA 1 Evaluation Panel") %>%
#   pivot_wider(names_from = assay_result_field, values_from = assay_result_value) %>%
#   filter(
#     # treatment_naive_at_visit == FALSE, 
#     ever_designated_as_elite_controller == FALSE#,
#     # days_since_eddi <= 1000,
#     # (eddi_interval_size <= 120 | days_since_eddi >= 2 * eddi_interval_size)
#   ) %>%
#   mutate(eddi = days_since_eddi, assay_value = final_result, viral_load = viral_load_closest_to_visit,
#          specimen = generic_result_identifier, subject_label_blinded = participant_identifier) %>%
#   group_by(subject_label_blinded) %>%
#   mutate(patients_ever_on_ART = max(treatment_naive_at_visit)) %>%
#   select(specimen, subject_label_blinded, viral_load, eddi, assay_value) %>%
#   select(subject_label_blinded, viral_load, eddi, ODn)


model1 <- nlme::lme(fixed = logvl ~ Group + bs(time_on_trt, 3),
    random = ~ 1 | subject_label_blinded,
    data = final_dataset,
    na.action = na.exclude, control = lmeControl(opt = "optim")
    )
summary(model1)
pred_vl <- predict(model1, newdata = final_dataset, level = 0:1)

model2 <- nlme::lme(fixed = logvl ~ Group + hiv_antibody_value + bs(time_on_trt, 3),
                    random = ~ 1 | subject_label_blinded,
                    data = final_dataset, # [final_dataset$time_on_trt >=0,]
                    na.action = na.exclude#, control = lmeControl(opt = "optim")
)
summary(model2)
pred_vl_2 <- predict(model2, newdata = final_dataset, level = 0:1)

model3 <- nlme::lme(fixed = logvl ~ hiv_antibody_value, #+ bs(time_on_trt, 3)
                    random = ~ 1 | subject_label_blinded,
                    data = final_dataset %>%
                      filter(time_on_trt>=0),
                    na.action = na.exclude, control = lmeControl(opt = "optim")
)
summary(model3)


data_without_sedia <- cbind((final_dataset %>%
  mutate(id = as.character(subject_label_blinded)) %>%
  select(id, subject_label_blinded, Group,time_on_trt, vl, logvl)), 
  pred_vl, pred_vl_2)