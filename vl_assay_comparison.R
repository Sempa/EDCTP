library(readxl)
library(tidyverse)
library(gtsummary)
library(dplyr)
library(ggplot2)
library(minpack.lm)
library(nlme)
library(splines)
library(cowplot)

source('estimate_measurement_noise_dist.R')
final_dataset <- read_csv("Sempa_final_pull_with_results.csv") %>%
  mutate(vl = ifelse(`Viral Load at Draw` == "<40", "40", ifelse(`Viral Load at Draw` == "NULL", "", `Viral Load at Draw`))) %>%
  mutate(
    logvl = log(as.numeric(vl), 10),
    time_on_trt = -1 * `days since tx start`,
    hiv_antibody_value = `Sedia LAg Odn screen`,
    RaceEthnicity = `Race/Ethnicity`
  ) %>%
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


x <- final_dataset[final_dataset$Group == "early suppressions", ] %>%
  select(
    Group, subject_label_blinded, `days since tx start`,
    `Viral Load at Draw`, vl, logvl, `Designated Elite at Draw`
  ) %>%
  arrange(subject_label_blinded, `days since tx start`)

model1 <- nlme::lme(
  fixed = logvl ~ Group + bs(time_on_trt, 3),
  random = ~ 1 | subject_label_blinded,
  data = final_dataset %>%
    filter(time_on_trt >= 0),
  na.action = na.exclude, control = lmeControl(opt = "optim")
)
summary(model1)

pred_vl <- predict(model1, newdata = final_dataset %>%
  filter(time_on_trt >= 0), level = 0:1)
after_ART <- final_dataset %>%
  mutate(subject_label_blinded = as.character(subject_label_blinded)) %>%
  filter(time_on_trt >= 0) %>%
  dplyr::select(subject_label_blinded, logvl, time_on_trt) %>%
  left_join(pred_vl %>%
    mutate(subject_label_blinded = as.character(subject_label_blinded)) %>%
    select(subject_label_blinded, predict.subject_label_blinded), by = "subject_label_blinded")

jpeg("website/EDCTP_website/vl_antibody_trajectory.jpeg", width = 480, height = 480, units = "px")
plot_grid(
  ggplot(
    data = final_dataset[final_dataset$Group == "early suppressions", ],
    aes(x = time_on_trt, y = logvl)
  ) +
    geom_line(aes(color = as.factor(subject_label_blinded), linetype = Group), size = 1.5) + # aes(color = subject_label_blinded, linetype = Group), group = Group
    geom_smooth(
      data = after_ART %>% filter(subject_label_blinded != "74498805"),
      aes(x = time_on_trt, y = predict.subject_label_blinded), size = 2,
      method = lm, formula = y ~ splines::bs(x, 3)
    ) +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = 18),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      legend.position = "none"
    ) +
    ylab(expression(paste(Log[10], " Viral load"))) +
    xlab("") +
    ggtitle("HIV viral load and HIV antibody response"),
  ggplot(
    data = final_dataset[final_dataset$Group == "early suppressions", ],
    aes(x = time_on_trt, y = `Sedia LAg Odn screen`)
  ) + # as.factor(subject_label_blinded)
    geom_line(aes(color = as.factor(subject_label_blinded), linetype = Group), size = 1.5) + # aes(color = subject_label_blinded, linetype = Group), group = Group
    geom_smooth(
      data = final_dataset[final_dataset$Group == "early suppressions", ] %>%
        filter(subject_label_blinded != 74498805) %>%
        filter(time_on_trt >= 0),
      aes(x = time_on_trt, y = `Sedia LAg Odn screen`), size = 2,
      method = lm, formula = y ~ splines::bs(x, 3)
    ) +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      legend.position = "none"
    ) +
    ylab("Sedia LAg ODn") +
    xlab("Time from ART start (Days)"),
  labels = "AUTO",
  ncol = 1,
  label_colour = "red"
)
# graphing the LAg results straddling ART start (the so-called `baseline` or last
# LAg reading before ART and the first LAg value at least 60months post ART start)
baseline_first_reading <- final_dataset %>%
  filter(time_on_trt <= 0) %>%
  group_by(subject_label_blinded) %>%
  mutate(last_LAg_reading_b4_ART = ifelse(time_on_trt <= 0, max(time_on_trt), 0), ) %>%
  filter(last_LAg_reading_b4_ART == time_on_trt) %>%
  dplyr::select(subject_label_blinded, time_on_trt) # , first_LAg_reading_on_ART, flag_visit

first_LAg_reading <- final_dataset %>%
  filter(time_on_trt > 180.5) %>%
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
  inner_join(visits_stradling_ART_start, by = c("subject_label_blinded", "time_on_trt"))

# graphing LAg reading straddling ART start

ggplot(
  data = LAg_stradling_ART_start,
  aes(x = time_on_trt, y = `Sedia LAg Odn screen`)
) + # as.factor(subject_label_blinded)
  geom_line(aes(color = as.factor(subject_label_blinded)), size = 1.5) + # aes(color = subject_label_blinded, linetype = Group), group = Group
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    legend.position = "none"
  ) +
  ylab("Sedia LAg ODn") +
  xlab("Time from ART start (Days)")

model_values_stradling_ARTstart <- nlme::lme(
  fixed = logvl ~ hiv_antibody_value,
  random = ~ 1 | subject_label_blinded,
  data = LAg_stradling_ART_start, # %>%
  # filter(time_on_trt>=0),
  na.action = na.exclude, control = lmeControl(opt = "optim")
)
summary(model_values_stradling_ARTstart)
#############################################################################################################
# studying the trajectory using CEPHIA other samples
sedia_generic <- read_csv("data/20180410-EP-LAgSedia-Generic.csv") %>%
  filter(visit_hivstatus == "P") %>%
  rename(sedia_ODn = result...72) %>%
  filter(visit_id != 21773 & visit_id != 21783 & visit_id != 21785)
data_generated <- sedia_generic %>%
  select(
    subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load, art_initiation_date, aids_diagnosis_date,
    art_interruption_date, art_resumption_date, treatment_naive,
    on_treatment, first_treatment
  ) %>%
  # arrange(subject_label_blinded, test_date)
  # filter(!is.na(art_initiation_date)) %>%
  arrange(subject_label_blinded, test_date) %>%
  # filter(!is.na(sedia_ODn)) %>%
  distinct(subject_label_blinded, test_date, .keep_all = T) %>%
  filter(!is.na(viral_load)) %>%
  group_by(subject_label_blinded) %>%
  # mutate(inter_test_interval = c(30, diff(test_date))) %>%
  # filter(inter_test_interval>=30) %>% #removing visits that are less than two months apart. CEPHIA has some visits that are even one day apart
  mutate(flagvl_100 = ifelse(viral_load<=100, 0,1),
         flagvl_400 = ifelse(viral_load<=400, 0,1),
         flagvl_1000 = ifelse(viral_load<=1000, 0,1)) %>%
  mutate(suprressed_throughout_followup_100 = ifelse(mean(flagvl_100) == 0 , 1, 0),
         suprressed_throughout_followup_400 = ifelse(mean(flagvl_400) ==0 , 1, 0),
         suprressed_throughout_followup_1000 = ifelse(mean(flagvl_1000) == 0 , 1, 0)
         ) %>% 
  # select(subject_label_blinded, test_date, viral_load,
  #        suprressed_throughout_followup_100, suprressed_throughout_followup_400, suprressed_throughout_followup_1000) %>%
  # filter(suprressed_throughout_followup_100==1)
  mutate(visits = 1:length(subject_label_blinded)) %>%
  mutate(n_visits = max(visits)) %>%
  filter(n_visits > 1) %>%
  ungroup() #%>%
  # mutate(
  #   time_on_trt = as.numeric(test_date - art_initiation_date),
  #   `Sedia LAg Odn screen` = sedia_ODn,
  #   vl_detectable_100 = (ifelse((viral_load) <= 100, 0, ifelse(viral_load > 100, 1, NA))),
  #   vl_detectable_400 = (ifelse((viral_load) <= 400, 0, ifelse(viral_load > 400, 1, NA))),
  #   vl_detectable_1000 = (ifelse((viral_load) <= 1000, 0, ifelse(viral_load > 1000, 1, NA)))
  # ) %>%
  # filter(time_on_trt >= 0) %>%
  # dplyr::select(
  #   subject_label_blinded, test_date, art_initiation_date, time_on_trt,
  #   `Sedia LAg Odn screen`, sedia_ODn, viral_load, n_visits, suprressed_throughout_followup, vl_detectable
  # ) %>% # , vl_detectable, inter_test_interval
  # arrange(subject_label_blinded, test_date)

ODn_regression_function <- function(dat) {
  # browser()
  model1 <- glm(
    formula = sedia_ODn_with_noise ~ days_since_eddi,
    family = gaussian(link = "identity"),
    data = dat
  )

  model2 <- glm(
    formula = sedia_ODn_with_noise ~ days_since_eddi,
    family = gaussian(link = "log"),
    data = dat
  )

  return(cbind(
    id = unique(dat$subject_label_blinded),
    intercept_link_identity = model1$coefficients[[1]],
    slope_link_identity = model1$coefficients[[2]],
    intercept_link_log = model2$coefficients[[1]],
    slope_link_log = model2$coefficients[[2]]
  ))
}


ggplot(
  data = data_generated %>%
    filter(suprressed_throughout_followup_100 == 1) %>%
    group_by(subject_label_blinded),
  aes(x = days_since_eddi, y = sedia_ODn)
) + # as.factor(subject_label_blinded)
  geom_line(aes(color = as.factor(subject_label_blinded)), size = 1) + # aes(color = subject_label_blinded, linetype = Group), group = Group
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    legend.position = "none"
  ) +
  ylab("Sedia LAg ODn") +
  xlab("Time since infection (Days)")

model_data_100 <- data.frame(id = NA, interceptlink_identity = NA, slope_link_identity = NA,
                         intercept_link_log = NA, slope_link_log = NA)
set.seed(11)
dat <- data_generated %>%
  filter(suprressed_throughout_followup_100 == 1) %>%
  mutate(sedia_ODn_with_noise = sedia_ODn + rnorm(n = length(subject_label_blinded), mean = 0, sd = noise$coefficients[1,1] + noise$coefficients[2,1] * sedia_ODn),
         sedia_ODn_with_noise = ifelse(sedia_ODn_with_noise<0,0, sedia_ODn_with_noise),
         sedia_ODn_with_noise = ifelse(is.na(sedia_ODn_with_noise), sedia_ODn, sedia_ODn_with_noise)
  ) %>%
  filter(!is.na(days_since_eddi)) %>%
  select(subject_label_blinded, test_date, days_since_eddi, sedia_ODn, sedia_ODn_with_noise, viral_load)
counter <- 0
for (i in 1:length(unique(dat$subject_label_blinded))) {
  counter <- counter + 1
  model_data_100[counter, ] <- ODn_regression_function(
    dat = subset(dat, dat$subject_label_blinded == unique(dat$subject_label_blinded)[i])
  )
  # print(i)
}
summary(model_data_100$slope_link_identity, na.rm = T)
model_data_100 %>%
  ggplot(aes(x = slope_link_identity)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, bins = 30, position = "identity") +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  labs(fill = "") +
  geom_density(
    lwd = 1.2,
    linetype = 1,
    colour = 2
  )

summary(model_data_100$slope_link_log, na.rm = T)
model_data_100 %>%
  ggplot(aes(x = slope_link_log)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, bins = 30, position = 'identity') +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  labs(fill = "") +
  geom_density(
    lwd = 1.2,
    linetype = 1,
    colour = 2
  )

################################################
#' Suppressed at 400
################################################

ggplot(
  data = data_generated %>%
    filter(suprressed_throughout_followup_400 == 1) %>%
    group_by(subject_label_blinded),
  aes(x = days_since_eddi, y = sedia_ODn)
) + # as.factor(subject_label_blinded)
  geom_line(aes(color = as.factor(subject_label_blinded)), size = 1) + # aes(color = subject_label_blinded, linetype = Group), group = Group
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    legend.position = "none"
  ) +
  ylab("Sedia LAg ODn") +
  xlab("Time since infection (Days)")

model_data_400 <- data.frame(id = NA, interceptlink_identity = NA, slope_link_identity = NA,
                             intercept_link_log = NA, slope_link_log = NA)
set.seed(11)
dat <- data_generated %>%
  filter(suprressed_throughout_followup_400 == 1) %>%
  mutate(sedia_ODn_with_noise = sedia_ODn + rnorm(n = length(subject_label_blinded), mean = 0, sd = noise$coefficients[1,1] + noise$coefficients[2,1] * sedia_ODn),
         sedia_ODn_with_noise = ifelse(sedia_ODn_with_noise<0,0, sedia_ODn_with_noise),
         sedia_ODn_with_noise = ifelse(is.na(sedia_ODn_with_noise), sedia_ODn, sedia_ODn_with_noise)
  ) %>%
  filter(!is.na(days_since_eddi)) %>%
  select(subject_label_blinded, test_date, days_since_eddi, sedia_ODn, sedia_ODn_with_noise, viral_load)
counter <- 0
for (i in 1:length(unique(dat$subject_label_blinded))) {
  counter <- counter + 1
  model_data_400[counter, ] <- ODn_regression_function(
    dat = subset(dat, dat$subject_label_blinded == unique(dat$subject_label_blinded)[i])
  )
  # print(i)
}
summary(model_data_400$slope_link_identity, na.rm = T)
model_data_400 %>%
  ggplot(aes(x = slope_link_identity)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, bins = 30, position = "identity") +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  labs(fill = "") +
  geom_density(
    lwd = 1.2,
    linetype = 1,
    colour = 2
  )

summary(model_data_400$slope_link_log, na.rm = T)
model_data_400 %>%
  ggplot(aes(x = slope_link_log)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, bins = 30, position = 'identity') +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  labs(fill = "") +
  geom_density(
    lwd = 1.2,
    linetype = 1,
    colour = 2
  )

################################################
#' Suppressed at 1000
################################################

ggplot(
  data = data_generated %>%
    filter(suprressed_throughout_followup_1000 == 1) %>%
    group_by(subject_label_blinded),
  aes(x = days_since_eddi, y = sedia_ODn)
) + # as.factor(subject_label_blinded)
  geom_line(aes(color = as.factor(subject_label_blinded)), size = 1) + # aes(color = subject_label_blinded, linetype = Group), group = Group
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    legend.position = "none"
  ) +
  ylab("Sedia LAg ODn") +
  xlab("Time since infection (Days)")

model_data_1000 <- data.frame(id = NA, interceptlink_identity = NA, slope_link_identity = NA,
                              intercept_link_log = NA, slope_link_log = NA)
set.seed(11)
dat <- data_generated %>%
  filter(suprressed_throughout_followup_1000 == 1) %>%
  mutate(sedia_ODn_with_noise = sedia_ODn + rnorm(n = length(subject_label_blinded), mean = 0, sd = noise$coefficients[1,1] + noise$coefficients[2,1] * sedia_ODn),
         sedia_ODn_with_noise = ifelse(sedia_ODn_with_noise<0,0, sedia_ODn_with_noise),
         sedia_ODn_with_noise = ifelse(is.na(sedia_ODn_with_noise), sedia_ODn, sedia_ODn_with_noise)
  ) %>%
  filter(!is.na(days_since_eddi)) %>%
  select(subject_label_blinded, test_date, days_since_eddi, sedia_ODn, sedia_ODn_with_noise, viral_load)
counter <- 0
for (i in 1:length(unique(dat$subject_label_blinded))) {
  counter <- counter + 1
  model_data_1000[counter, ] <- ODn_regression_function(
    dat = subset(dat, dat$subject_label_blinded == unique(dat$subject_label_blinded)[i])
  )
  # print(i)
}
summary(model_data_1000$slope_link_identity, na.rm = T)
model_data_1000 %>%
  ggplot(aes(x = slope_link_identity)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, bins = 30, position = "identity") +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  labs(fill = "") +
  geom_density(
    lwd = 1.2,
    linetype = 1,
    colour = 2
  )

summary(model_data_1000$slope_link_log, na.rm = T)
model_data_1000 %>%
  ggplot(aes(x = slope_link_log)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, bins = 30, position = 'identity') +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  labs(fill = "") +
  geom_density(
    lwd = 1.2,
    linetype = 1,
    colour = 2
  )
