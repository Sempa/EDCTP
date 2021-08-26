library(readxl)
library(tidyverse)
sheet1 <- read_excel("Sempa_query_withLAgresults.xlsx", 
                     sheet = "early_suppressions") %>%
  arrange(subject_label_blinded, `days since tx start`)
sheet2 <- read_excel("Sempa_query_withLAgresults.xlsx", 
                     sheet = "suppression_failures")%>%
  arrange(subject_label_blinded, `days since tx start`)
Sempa_query_withLAgresults <- rbind(sheet1, sheet2)
sedia_missing <-Sempa_query_withLAgresults %>%
  filter(is.na(`LAg result?`) | viral_load=='NULL')
write.csv(sedia_missing, 'specimens_to_be_tested.csv')
