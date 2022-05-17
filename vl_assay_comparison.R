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

model1 <- nlme::lme(fixed = logvl ~ Group + bs(time_on_trt, 3),
                    random = ~ 1 | subject_label_blinded,
                    data = final_dataset %>%
                      filter(time_on_trt>=0),
                    na.action = na.exclude, control = lmeControl(opt = "optim")
)
summary(model1)

pred_vl <- predict(model1, newdata = final_dataset %>%
                     filter(time_on_trt>=0), level = 0:1)
after_ART <- final_dataset %>%
  mutate(subject_label_blinded = as.character(subject_label_blinded)) %>%
  filter(time_on_trt>=0) %>%
  dplyr::select(subject_label_blinded, logvl, time_on_trt) %>%
  left_join(pred_vl %>% 
              mutate(subject_label_blinded = as.character(subject_label_blinded)) %>% 
              select(subject_label_blinded, predict.subject_label_blinded), by = 'subject_label_blinded')

jpeg('website/EDCTP_website/vl_antibody_trajectory.jpeg', width = 480, height = 480, units = "px")
plot_grid(
ggplot(data = final_dataset[final_dataset$Group == 'early suppressions',],
       aes(x = time_on_trt, y = logvl)) +
  geom_line(aes(color = as.factor(subject_label_blinded), linetype = Group), size = 1.5) +#aes(color = subject_label_blinded, linetype = Group), group = Group
  geom_smooth(data = after_ART %>% filter(subject_label_blinded!='74498805'), 
              aes(x = time_on_trt, y = predict.subject_label_blinded), size = 2,
              method = lm, formula = y ~ splines::bs(x, 3)) +
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
  ylab(expression(paste(Log[10], ' Viral load'))) + xlab('') +
  ggtitle('HIV viral load and HIV antibody response'),

ggplot(data = final_dataset[final_dataset$Group == 'early suppressions',],
       aes(x = time_on_trt, y = `Sedia LAg Odn screen`)) + #as.factor(subject_label_blinded)
  geom_line(aes(color = as.factor(subject_label_blinded), linetype = Group), size = 1.5) +#aes(color = subject_label_blinded, linetype = Group), group = Group
  geom_smooth(data = final_dataset[final_dataset$Group == 'early suppressions',] %>% 
                filter(subject_label_blinded!=74498805) %>%
                filter(time_on_trt>=0), 
              aes(x = time_on_trt, y = `Sedia LAg Odn screen`), size = 2,
              method = lm, formula = y ~ splines::bs(x, 3)) +
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

model_values_stradling_ARTstart <- nlme::lme(fixed = logvl ~ hiv_antibody_value,
                    random = ~ 1 | subject_label_blinded,
                    data = LAg_stradling_ART_start, #%>%
                      # filter(time_on_trt>=0),
                    na.action = na.exclude, control = lmeControl(opt = "optim")
)
summary(model_values_stradling_ARTstart)
#############################################################################################################
# studying the trajectory using CEPHIA samples that don't have pre 
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
  filter(!is.na(viral_load)) %>%
  group_by(subject_label_blinded) %>%
  # mutate(inter_test_interval = c(30, diff(test_date))) %>%
  # filter(inter_test_interval>=30) %>% #removing visits that are less than two months apart. CEPHIA has some visits that are even one day apart
  mutate(mean_vl = mean(viral_load)) %>%
  mutate(suprressed_throughout_followup = ifelse(mean_vl<=40, 1,0)) %>%
  mutate(visits = 1:length(subject_label_blinded)) %>%
  mutate(n_visits = max(visits)) %>%
  filter(n_visits>1) %>%
  ungroup() %>%
  mutate(time_on_trt = as.numeric(test_date - art_initiation_date),
         `Sedia LAg Odn screen` = sedia_ODn,
         vl_detectable = (ifelse((viral_load) <=40, 0, ifelse(viral_load >40, 1, NA)))
         ) %>%
  dplyr::select(subject_label_blinded, test_date, art_initiation_date, time_on_trt, 
                `Sedia LAg Odn screen`, sedia_ODn, viral_load, n_visits, suprressed_throughout_followup, vl_detectable) %>% #, vl_detectable, inter_test_interval
  arrange(subject_label_blinded, test_date)

model_cephia_2 <- nlme::lme(fixed = log10(viral_load) ~ sedia_ODn, #+ bs(time_on_trt, 3)
                           random = ~ 1 | subject_label_blinded,
                           data = data_generated %>%
                             filter(time_on_trt>=0),
                           na.action = na.exclude, control = lmeControl(opt = "optim")
)
summary(model_cephia_2)

linear_reg_function <- function(dat, type) {
  browser()
  if (type == 1) {
    if (dat$suprressed_throughout_followup == 0) {
      model <- glm(
        formula = log10(viral_load) ~ sedia_ODn,
        family = gaussian,
        data = dat
      )
    } else if (dat$suprressed_throughout_followup == 1) {
      model <- glm(
        formula = sedia_ODn ~ time_on_trt,
        family = gaussian,
        data = dat
      )
    }
    return(cbind(
      id = unique(dat$subject_label_blinded),
      intercept = model$coefficients[[1]],
      slope = model$coefficients[[2]],
      suppressed = unique(dat$suprressed_throughout_followup)
    ))
  } else if (type == 2) {
    if (dat$vl_detectable == 1) {
      model <- glm(
        formula = log10(viral_load) ~ sedia_ODn,
        family = gaussian,
        data = dat
      )
    } else if (dat$vl_detectable == 0) {
      model <- glm(
        formula = sedia_ODn ~ time_on_trt,
        family = gaussian,
        data = dat
      )
    }
    return(cbind(
      id = unique(dat$subject_label_blinded),
      intercept = model$coefficients[[1]],
      slope = model$coefficients[[2]],
      suppressed = unique(dat$vl_detectable)
    ))
  }
}
#####dealing with unsuppressed visits
model_data <- data.frame(id=NA, intercept = NA, slope = NA, suppressed=NA)
counter <- 0
unsuppressed_visit <- data_generated %>%
  filter(vl_detectable==1)
for (i in 1:length(unique(unsuppressed_visit$subject_label_blinded))) {
  counter <- counter + 1
  model_data[counter,] <- linear_reg_function(dat = subset(unsuppressed_visit, unsuppressed_visit$subject_label_blinded== unique(unsuppressed_visit$subject_label_blinded)[i]),
                                              type = 2)
}
mean(model_data$slope, na.rm = T)
summary(lm(slope~1, data = model_data))

#####dealing with suppressed visits
model_data_supp <- data.frame(id=NA, intercept = NA, slope = NA, suppressed=NA)
counter <- 0
suppressed_visit <- data_generated %>%
  filter(vl_detectable==0)
for (i in 1:length(unique(suppressed_visit$subject_label_blinded))) {
  counter <- counter + 1
  model_data_supp[counter,] <- linear_reg_function(dat = subset(suppressed_visit, suppressed_visit$subject_label_blinded== unique(suppressed_visit$subject_label_blinded)[i]),
                                              type = 2)
}

mean(model_data_supp$slope, na.rm = T)
summary(lm(slope~1, data = model_data_supp))

model_data_supp%>%
  ggplot(aes(x=slope, type=suppressed)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="")
cols <- c("#F76D5E", "#FFFFBF")#, "#72D8FF"
v=rbind(model_data, model_data_supp) %>%
  filter(!is.na(slope)) %>%
  # filter(slope>-10) %>%
  # filter(slope<5) %>%
  ggplot(aes(x = slope, colour = as.factor(suppressed))) +
    geom_density(lwd = 1.2, linetype = 1) + 
  scale_color_manual(values = cols)
v


after_ART <- final_dataset %>%
  mutate(subject_label_blinded = as.character(subject_label_blinded)) %>%
  filter(time_on_trt>=0) %>%
  dplyr::select(subject_label_blinded, logvl, time_on_trt) %>%
  left_join(pred_vl %>% 
              mutate(subject_label_blinded = as.character(subject_label_blinded)) %>% 
              select(subject_label_blinded, predict.subject_label_blinded), by = 'subject_label_blinded') %>%
  mutate(predict_exp_fit = exp(lambda$parameters[[1]]) * time_on_trt)

ggplot(data = final_dataset[final_dataset$Group == 'early suppressions',],
       aes(x = time_on_trt, y = logvl)) +
  geom_line(aes(color = as.factor(subject_label_blinded), linetype = Group), size = 1.5) +#aes(color = subject_label_blinded, linetype = Group), group = Group
  geom_smooth(data = after_ART %>% filter(subject_label_blinded!='74498805'), 
              aes(x = time_on_trt, y = predict.subject_label_blinded), size = 2,
              method = lm, formula = y ~ splines::bs(x, 3)) +
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
  ylab(expression(paste(Log[10], ' Viral load'))) + xlab('')


model2 <- nlme::lme(fixed = logvl ~ Group + hiv_antibody_value + bs(time_on_trt, 3),
                    random = ~ 1 | subject_label_blinded,
                    data = final_dataset%>%
                      filter(time_on_trt>=0), # [final_dataset$time_on_trt >=0,]
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