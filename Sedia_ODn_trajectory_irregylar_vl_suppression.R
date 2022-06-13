
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
  mutate(flagvl = ifelse(viral_load > 1000, 0, 1)) %>%
  mutate(irregular_suprression = ifelse(mean(flagvl) == 0, 1, 0)) %>%
  filter(irregular_suprression == 0) %>%
  mutate(visits = 1:length(subject_label_blinded)) %>%
  mutate(n_visits = max(visits)) %>%
  filter(n_visits > 1) %>%
  arrange(subject_label_blinded, days_since_eddi) %>%
  ungroup()


#' does the slope for ODn among those with ODn<2 differ from those with ODn>2 in
#' patients peaking?
write.csv(data_intermitent_suppression, "output_table/intermitent_suppression.csv")
data_intermitent_suppression_selected_visits <- read.csv("output_table/intermitent_suppression_selected.csv")
data_intermitent_suppression_toPeak_visits <- data_intermitent_suppression_selected_visits %>%
  filter(!is.na(to_peak)) %>%
  unite("id", c(subject_label_blinded, set_to_peak))

data_intermitent_suppression_totrough_visits <- data_intermitent_suppression_selected_visits %>%
  filter(!is.na(to_trough)) %>%
  unite("id", c(subject_label_blinded, set_to_trough))
summary(to_peak <- data_intermitent_suppression_selected_visits %>%
  filter(!is.na(to_peak)) %>%
  unite("id", c(subject_label_blinded, set_to_peak)) %>%
  group_by(id) %>%
  summarise(x = mean(diff(days_since_eddi))))

summary(to_trough <- data_intermitent_suppression_selected_visits %>%
  filter(!is.na(to_trough)) %>%
  unite("id", c(subject_label_blinded, set_to_trough)) %>%
  group_by(id) %>%
  summarise(x = mean(diff(days_since_eddi))))

slopes_for_unsuppressed <- function(ODn_vl_data, threshold) {
  # browser()
  dir.create(paste("unsupp_", threshold, Sys.Date(), sep = "_"))
  folder_name <- paste("unsupp_", threshold, Sys.Date(), sep = "_")
  # data_generated <- ODn_vl_data %>%
  #   select(
  #     subject_label_blinded, days_since_eddi, test_date, sedia_ODn, viral_load # , art_initiation_date, aids_diagnosis_date,
  #     # art_interruption_date, art_resumption_date, treatment_naive,
  #     # on_treatment, first_treatment
  #   ) %>%
  #   arrange(subject_label_blinded, test_date) %>%
  #   # distinct(subject_label_blinded, test_date, .keep_all = T) %>%
  #   filter(!is.na(viral_load)) %>%
  #   group_by(subject_label_blinded) %>%
  #   mutate(flagvl = ifelse(viral_load <= threshold, 0, 1)) %>%
  #   mutate(suprressed_throughout_followup = ifelse(mean(flagvl) == 0, 1, 0)) %>%
  #   mutate(visits = 1:length(subject_label_blinded)) %>%
  #   mutate(n_visits = max(visits)) %>%
  #   filter(n_visits > 1) %>%
  #   ungroup()

  jpeg(paste(folder_name, "/", "Figure1", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)

  print(
    ggplot(
      data = ODn_vl_data %>%
        # filter(suprressed_throughout_followup == 1) %>%
        group_by(id),
      aes(x = days_since_eddi, y = sedia_ODn)
    ) + # as.factor(id)
      geom_line(aes(color = as.factor(id)), size = 1) + # aes(color = id, linetype = Group), group = Group
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
      xlab("Time since infection (Days)") # +
    # ggtitle(paste("Sedia ODn vs EDDI_", threshold))
  )
  dev.off()

  model_data <- data.frame(
    id = NA, interceptlink_identity = NA, slope_link_identity = NA,
    intercept_link_log = NA, slope_link_log = NA
  )
  # set.seed(11)
  dat <- ODn_vl_data %>%
    # filter(suprressed_throughout_followup == 1) %>%
    # mutate(
    #   sedia_ODn_with_noise = sedia_ODn + rnorm(n = length(id), mean = 0, sd = noise$coefficients[1, 1] + noise$coefficients[2, 1] * sedia_ODn),
    #   sedia_ODn_with_noise = ifelse(sedia_ODn_with_noise < 0, 0, sedia_ODn_with_noise),
    #   sedia_ODn_with_noise = ifelse(is.na(sedia_ODn_with_noise), sedia_ODn, sedia_ODn_with_noise)
    # ) %>%
    filter(!is.na(days_since_eddi)) # %>%
  # select(id, test_date, days_since_eddi, sedia_ODn, sedia_ODn_with_noise, viral_load)
  counter <- 0
  for (i in 1:length(unique(dat$id))) {
    counter <- counter + 1
    model_data[counter, ] <- ODn_regression_function(
      dat = subset(dat, dat$id == unique(dat$id)[i])
    )
    # print(i)
  }
  # browser()
  # summary(model_data$slope_link_identity, na.rm = T)

  jpeg(paste(folder_name, "/", "Figure_link_identity", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)

  print(
    model_data %>%
      mutate(slope_link_identity = as.numeric(slope_link_identity)) %>%
      ggplot(aes(x = slope_link_identity)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, bins = 30, position = "identity") +
      scale_fill_manual(values = c("#69b3a2", "#404080")) +
      labs(fill = "") +
      geom_density(
        lwd = 1.2,
        linetype = 1,
        colour = 2
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
      ggtitle(paste("histogram mean slope link=identity_", threshold))
  )
  dev.off()
  # summary(model_data$slope_link_log, na.rm = T)

  jpeg(paste(folder_name, "/", "Figure_link_log", threshold, ".jpeg", sep = ""), units = "in", width = 8, height = 6, res = 300)

  print(
    model_data %>%
      mutate(slope_link_log = as.numeric(slope_link_log)) %>%
      ggplot(aes(x = slope_link_log)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, bins = 30, position = "identity") +
      scale_fill_manual(values = c("#69b3a2", "#404080")) +
      labs(fill = "") +
      geom_density(
        lwd = 1.2,
        linetype = 1,
        colour = 2
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
      ggtitle(paste("histogram mean slope link=log_", threshold))
  )
  dev.off()

  return(
    cbind(
      suppression = threshold,
      slope_mean_link_identity = mean(as.numeric(model_data$slope_link_identity)),
      pvalue_link_identity_slope_diff_zero = t.test((as.numeric(model_data$slope_link_identity)))$p.value,
      slope_mean_link_log = mean(as.numeric(model_data$slope_link_log)),
      pvalue_link_log_slope_diff_zero = t.test((as.numeric(model_data$slope_link_log)))$p.value
    )
  )
}

ODn_regression_function <- function(dat) {
  # browser()
  model1 <- glm(
    formula = sedia_ODn ~ days_since_eddi,
    family = gaussian(link = "identity"),
    data = dat
  )

  model2 <- glm(
    formula = sedia_ODn ~ days_since_eddi,
    family = gaussian(link = "log"),
    data = dat
  )

  return(cbind(
    id = unique(dat$id),
    intercept_link_identity = model1$coefficients[[1]],
    slope_link_identity = model1$coefficients[[2]],
    intercept_link_log = model2$coefficients[[1]],
    slope_link_log = model2$coefficients[[2]]
  ))
}


x <- slopes_for_unsuppressed(ODn_vl_data = data_intermitent_suppression_toPeak_visits, threshold = "to_unsuppressed")
x1 <- slopes_for_unsuppressed(ODn_vl_data = data_intermitent_suppression_totrough_visits, threshold = "from_unsuppressed")
unsuppressed_group_slopes <- data.frame(rbind(x, x1))