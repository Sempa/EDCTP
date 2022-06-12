data_intermitent_suppression <- sedia_generic %>%
  select(
    subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load # , art_initiation_date, aids_diagnosis_date,
    # art_interruption_date, art_resumption_date, treatment_naive,
    # on_treatment, first_treatment
  ) %>%
  arrange(subject_label_blinded, test_date) %>%
  # distinct(subject_label_blinded, test_date, .keep_all = T) %>%
  filter(!is.na(viral_load)) %>%
  group_by(subject_label_blinded) %>%
  mutate(flagvl = ifelse(viral_load >1000, 0, 1)) %>%
  mutate(irregular_suprression = ifelse(mean(flagvl) == 0, 1, 0)) %>%
  filter(irregular_suprression == 0) %>%
  mutate(visits = 1:length(subject_label_blinded)) %>%
  mutate(n_visits = max(visits)) %>%
  filter(n_visits > 1) %>%
  arrange(subject_label_blinded, test_date) %>%
  ungroup()


#' does the slope for ODn among those with ODn<2 differ from those with ODn>2 in 
#' patients peaking? 
write.csv(data_intermitent_suppression, 'output_table/intermitent_suppression.csv')
