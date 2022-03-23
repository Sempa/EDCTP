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

dev.off()

# # after_ART <- final_dataset %>%
# #   mutate(subject_label_blinded = as.character(subject_label_blinded)) %>%
# #   filter(time_on_trt>=0) %>%
# #   dplyr::select(subject_label_blinded, logvl, time_on_trt) %>%
# #   left_join(pred_vl %>% 
# #               mutate(subject_label_blinded = as.character(subject_label_blinded)) %>% 
# #               select(subject_label_blinded, predict.subject_label_blinded), by = 'subject_label_blinded')
# ## Exponential model
# library(minpack.lm)
# nls.fit <- function(lambda) {
#   modle1 <- nlsLM(logvl ~ exp(-lambda * time_on_trt),
#                   data = final_dataset %>%
#                     filter(time_on_trt>=0),
#                   start = list(lambda = lambda)
#   ) # , verbose = FALSE)
#   return(modle1)
# }
# lambda <- summary(nls.fit(0))

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