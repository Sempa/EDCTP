library(tidyverse)
library(readr)
library(splines)
pt_data <- read_csv("Sempa_final_pull_with_results.csv") %>%
  mutate(vl = ifelse(`Viral Load at Draw`== '<40', '40', ifelse(`Viral Load at Draw`=='NULL', "", `Viral Load at Draw`))) %>%
  mutate(logvl = log(as.numeric(vl),10),
         time_on_trt = -1 * `days since tx start`,
         hiv_antibody_value = `Sedia LAg Odn screen`) %>%
  arrange(subject_label_blinded, time_on_trt)

dt01 <- pt_data %>%
  distinct(subject_label_blinded, .keep_all = TRUE) %>%
  group_by(BiologicalSex) %>%
  summarise(`Freq.` = n()) %>%
  mutate(
    percentage = round((`Freq.` / (sum(`Freq.`))) * 100, 1)
)
dt02 <- pt_data %>%
  distinct(subject_label_blinded, .keep_all = TRUE) %>%
  group_by(Subtype_confirmed) %>%
  summarise(`Freq.` = n()) %>%
  mutate(
    percentage = round((`Freq.` / (sum(`Freq.`))) * 100, 1)
  )

dt03 <- pt_data %>%
  distinct(subject_label_blinded, .keep_all = TRUE) %>%
  summarise(
    N = n(), Average = round(mean(`days since tx start`), 2), 
    `standard devaition` = round(sd(`days since tx start`), 2), 
    `median value` = median(`days since tx start`),
    `lower quartile` = quantile(`days since tx start`, .25), 
    `upper quartile` = quantile(`days since tx start`, .75)
  )

x=pt_data[pt_data$Group == 'early suppressions',] %>% 
  select(Group, subject_label_blinded, `days since tx start`, 
         `Viral Load at Draw`, vl, logvl, `Designated Elite at Draw`) %>% 
  arrange(subject_label_blinded, `days since tx start`)

ggplot(data = pt_data[pt_data$Group == 'early suppressions',],
       aes(x = time_on_trt, y = logvl)) + #as.factor(subject_label_blinded)
  geom_line(aes(color = as.factor(subject_label_blinded), linetype = Group), size = 1.5) +#aes(color = subject_label_blinded, linetype = Group), group = Group
  theme_bw() +
  theme(legend.position = "none")

ggplot(data = pt_data[pt_data$Group == 'early suppressions',],
       aes(x = time_on_trt, y = `Sedia LAg Odn screen`)) + #as.factor(subject_label_blinded)
  geom_line(aes(color = as.factor(subject_label_blinded), linetype = Group), size = 1.5) +#aes(color = subject_label_blinded, linetype = Group), group = Group
  theme_bw() +
  theme(legend.position = "none")
  
model1 <- nlme::lme(fixed = logvl ~ Group + bs(time_on_trt, 3),
    random = ~ 1 | subject_label_blinded,
    data = pt_data,
    na.action = na.exclude#, control = lmeControl(opt = "optim")
    )
summary(model1)
pred_vl <- predict(model1, newdata = pt_data, level = 0:1)

model2 <- nlme::lme(fixed = logvl ~ Group + hiv_antibody_value + bs(time_on_trt, 3),
                    random = ~ 1 | subject_label_blinded,
                    data = pt_data, # [pt_data$time_on_trt >=0,]
                    na.action = na.exclude#, control = lmeControl(opt = "optim")
)
summary(model2)
pred_vl_2 <- predict(model2, newdata = pt_data, level = 0:1)

data_without_sedia <- cbind((pt_data %>%
  mutate(id = as.character(subject_label_blinded)) %>%
  select(id, subject_label_blinded, Group,time_on_trt, vl, logvl)), 
  pred_vl, pred_vl_2)