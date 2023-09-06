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
         `ODs 3` = `...14`) %>%
  dplyr::select(id, `Date Collected`, `Study Number`, `ODs 1`, `ODs 2`,`ODs 3`) %>%
  arrange(`Study Number`, `Date Collected`)
dt01 <- africos_data %>%
  pivot_longer(cols = c('ODs 1', 'ODs 2','ODs 3'),
               names_to = 'ODs',
               values_to = 'OD')
# length(unique(dt01$`Study Number`))
