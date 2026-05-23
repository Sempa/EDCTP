library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(plotly)
library(tidyverse)
library(htmlwidgets)
library(webshot2)

# --- Parameter grid ---
spec_values <- c(seq(0.60, 0.90, 0.10), 0.95)
AB_rebound_delay_seq <- seq(0, 0.4, by = 0.01)
annual_rebound_rate_seq <- seq(0.04, 0.5, by = 0.01)

param_grid <- expand.grid(
  AB_Specificity = spec_values,
  AB_rebound_delay = AB_rebound_delay_seq,
  annual_rebound_rate = annual_rebound_rate_seq
)

population <- 1000

# --- Function that computes metrics for one parameter set ---
compute_metrics <- function(AB_Specificity, AB_rebound_delay, annual_rebound_rate, population = 1000) {
  monitoring_scenario <- c('annual VL', 'Annual AB', 'biannual AB')
  N <- rep(population, length(monitoring_scenario))
  prim_VL_tests <- c(population, 0, 0)
  testing_interval <- c(1, 1, 0.5)
  AB_testing <- c(0, rep(population, 2))
  
  number_of_rebounds_PerAnnum <- rep(population * annual_rebound_rate, 3)
  number_of_rebounds_Per_test_round <- c(
    NA,
    population * annual_rebound_rate * (testing_interval[2] + AB_rebound_delay),
    (testing_interval[3] * annual_rebound_rate * population) + (population * annual_rebound_rate * AB_rebound_delay)
  )
  
  true_positives <- c(
    population * annual_rebound_rate,
    number_of_rebounds_PerAnnum[2] * testing_interval[2],
    number_of_rebounds_PerAnnum[3] * testing_interval[3]
  )
  
  false_positives <- c(
    0,
    (N[2] - number_of_rebounds_PerAnnum[2]) * (1 - AB_Specificity),
    (N[3] - number_of_rebounds_PerAnnum[3]) * (1 - AB_Specificity)
  )
  
  vl_confirmation <- c(
    0,
    true_positives[2] + false_positives[2],
    true_positives[3] + false_positives[3]
  )
  
  total_VL <- c(
    vl_confirmation[1] + prim_VL_tests[1],
    vl_confirmation[2] + prim_VL_tests[2],
    vl_confirmation[3] + prim_VL_tests[3]
  )
  
  mean_delay_yrs <- c(
    0.5,
    0.5 * (AB_rebound_delay + testing_interval[2] + AB_rebound_delay),
    0.5 * (AB_rebound_delay + testing_interval[3] + AB_rebound_delay)
  )
  mean_delay_months <- mean_delay_yrs * 12
  
  total_VL_per_yr <- c(
    N[1],
    total_VL[2] / testing_interval[2],
    total_VL[3] / testing_interval[3]
  )
  
  tibble(
    Scenario = monitoring_scenario,
    AB_Specificity = AB_Specificity,
    AB_rebound_delay = AB_rebound_delay,
    annual_rebound_rate = annual_rebound_rate,
    Mean_delay_months = mean_delay_months,
    Total_VL_per_year = total_VL_per_yr
  )
}

# --- Run across grid ---
results <- param_grid %>%
  mutate(data = pmap(
    list(AB_Specificity, AB_rebound_delay, annual_rebound_rate),
    ~ compute_metrics(..1, ..2, ..3, population)
  )) %>%
  unnest(data, names_sep = "_") %>%
  mutate(Total_VL_saved = (population - data_Total_VL_per_year)/population) %>%
  dplyr::select(
    Scenario = data_Scenario, 
    AB_Specificity, AB_rebound_delay, annual_rebound_rate, 
    Mean_delay_months = data_Mean_delay_months, 
    Total_VL_per_year = data_Total_VL_per_year,
    Total_VL_saved
  )

# --- Plot distributions ---

x <- results %>%
  filter(Scenario == 'Annual AB')

p1 <- plotly::plot_ly(
  data = x,
  x = ~AB_Specificity,
  y = ~annual_rebound_rate,
  z = ~Total_VL_saved,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 4, color = ~Total_VL_saved, colorscale = "Viridis")
) %>%
  layout(
    title = list(
      # text = "Percentage total savings on VL testing: Annual AB testing",
      x = 0.5,                # center title like ggplot hjust = 0.5
      font = list(size = 22)  # title font size
    ),
    scene = list(
      xaxis = list(
        title = list(text = "AB Specificity", font = list(size = 18)),
        tickfont = list(size = 16),
        showbackground = FALSE,  # match panel.background = element_blank()
        showgrid = TRUE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = list(text = "Annual Rebound Rate", font = list(size = 18)),
        tickfont = list(size = 16),
        showbackground = FALSE,
        showgrid = TRUE,
        zeroline = FALSE
      ),
      zaxis = list(
        title = list(text = "Total VL Saved", font = list(size = 18)),
        tickfont = list(size = 16),
        range = c(0, max(x$Total_VL_saved, na.rm = TRUE)),
        showbackground = FALSE,
        showgrid = TRUE,
        zeroline = FALSE
      ),
      aspectmode = "manual",
      aspectratio = list(x = 1, y = 1, z = 0.6),
      camera = list(eye = list(x = 1.6, y = 1.6, z = 0.9))
    ),
    margin = list(l = 10, r = 10, b = 10, t = 60),
    paper_bgcolor = "white",
    plot_bgcolor = "white",
    font = list(family = "sans-serif", size = 18) # overall text size
  ) %>%
  plotly::layout(
    # title = list(text = "3D Points: Total VL Saved", x = 0.5, font = list(size = 36)),
    font = list(size = 20),     # base text size
    scene = list(
      xaxis = list(title = list(text = "AB Specificity", font = list(size = 20)), tickfont = list(size = 18)),
      yaxis = list(title = list(text = "Annual Rebound Rate", font = list(size = 20)), tickfont = list(size = 18)),
      zaxis = list(title = list(text = "Total VL Saved", font = list(size = 20)), tickfont = list(size = 18))
    )
  )
p1
# plotly::save_image(p1, "annual AB.png", width = 1600, height = 1200, scale = 2)
# htmlwidgets::saveWidget(p1, "biannual_AB.html", selfcontained = TRUE)
# webshot2::webshot("biannual_AB.html", "biannual_AB.png", 
#                   vwidth = 1800, vheight = 1400, zoom = 2)
x <- results %>%
  filter(Scenario == 'biannual AB')

p2 <- plotly::plot_ly(
  data = x,
  x = ~AB_Specificity,
  y = ~annual_rebound_rate,
  z = ~Total_VL_saved,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 4, color = ~Total_VL_saved, colorscale = "Viridis")
) %>%
  layout(
    title = list(
      # text = "Percentage total savings on VL testing: Biannual AB testing",
      x = 0.5,                # center title like ggplot hjust = 0.5
      font = list(size = 22)  # title font size
    ),
    scene = list(
      xaxis = list(
        title = list(text = "AB Specificity", font = list(size = 18)),
        tickfont = list(size = 18),
        showbackground = FALSE,  # match panel.background = element_blank()
        showgrid = TRUE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = list(text = "Annual Rebound Rate", font = list(size = 18)),
        tickfont = list(size = 18),
        showbackground = FALSE,
        showgrid = TRUE,
        zeroline = FALSE
      ),
      zaxis = list(
        title = list(text = "Total VL Saved", font = list(size = 18)),
        tickfont = list(size = 18),
        tickvals = c(0, 0.2, 0.4, 0.6, 0.8, 1),
        ticktext = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"),
        range = c(0, max(x$Total_VL_saved, na.rm = TRUE)),
        showbackground = FALSE,
        showgrid = TRUE,
        zeroline = FALSE
      ),
      aspectmode = "manual",
      aspectratio = list(x = 1, y = 1, z = 0.6),
      camera = list(eye = list(x = 1.6, y = 1.6, z = 0.9))
    ),
    margin = list(l = 10, r = 10, b = 10, t = 60),
    paper_bgcolor = "white",
    plot_bgcolor = "white",
    font = list(family = "sans-serif", size = 20) # overall text size
  ) %>%
  plotly::layout(
    # title = list(text = "3D Points: Total VL Saved", x = 0.5, font = list(size = 36)),
    font = list(size = 20),     # base text size
    scene = list(
      xaxis = list(title = list(text = "AB Specificity", font = list(size = 20)), tickfont = list(size = 18)),
      yaxis = list(title = list(text = "Annual Rebound Rate", font = list(size = 20)), tickfont = list(size = 18)),
      zaxis = list(title = list(text = "Total VL Saved", font = list(size = 20)), tickfont = list(size = 18))
    )
  )
p2
# plotly::save_image(p2, "biannual AB.png", width = 1600, height = 1200, scale = 2)
htmlwidgets::saveWidget(p2, "biannual_AB.html", selfcontained = TRUE)
webshot2::webshot("biannual_AB.html", "biannual_AB.png", 
                  vwidth = 1800, vheight = 1400, zoom = 2)
## 1. Distribution of Mean Delay (months)
p3 <- ggplot(results %>%
               filter(Scenario == 'Annual AB' & AB_Specificity == 0.6 & annual_rebound_rate == 0.04), 
             aes(y = Mean_delay_months, x = AB_rebound_delay)) + #, fill = Scenario
  geom_point(alpha = 0.6, size = 5) +
  labs(
    # title = "Distribution of Mean Delay (months)",
    # subtitle = "Varying AB rebound delay (0–0.4) and annual rebound rate (0.04–0.5)",
    x = "AB rebound delay",
    y = "Mean Delay (months)"#,
    # fill = "Scenario"
  ) +
  theme(
    text = element_text(size = 22),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "null")
  ) +
  scale_fill_brewer(palette = "Set2")

ggsave("epidemics/mean_delay_annual_AB.png", plot = p3,
       width = 10, height = 6, dpi = 300)

## 2. Distribution of Total VL (per year)

p4 <- ggplot(results %>%
               filter(Scenario == 'biannual AB' & AB_Specificity == 0.6 & annual_rebound_rate == 0.04), 
             aes(y = Mean_delay_months, x = AB_rebound_delay)) + #, fill = Scenario
  geom_point(alpha = 0.6, size = 5) +
  labs(
    # title = "Distribution of Mean Delay (months)",
    # subtitle = "Varying AB rebound delay (0–0.4) and annual rebound rate (0.04–0.5)",
    x = "AB rebound delay",
    y = "Mean Delay (months)"#,
    # fill = "Scenario"
  ) +
  theme(
    text = element_text(size = 22),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "null")
  ) +
  scale_fill_brewer(palette = "Set2")

ggsave("epidemics/mean_delay_biannual_AB.png", plot = p4,
       width = 10, height = 6, dpi = 300)


# -----------------------------
# Representative scenarios
# -----------------------------
rep_specs <- c(0.80, 0.90, 0.95)
rep_delays <- c(0.00, 0.10, 0.20)
rep_rebounds <- c(0.05, 0.10, 0.20)

# Because of floating point issues, round first
results_rep <- results %>%
  mutate(
    AB_Specificity = round(AB_Specificity, 2),
    AB_rebound_delay = round(AB_rebound_delay, 2),
    annual_rebound_rate = round(annual_rebound_rate, 2)
  ) %>%
  filter(
    Scenario %in% c("Annual AB", "biannual AB"),
    AB_Specificity %in% rep_specs,
    AB_rebound_delay %in% rep_delays,
    annual_rebound_rate %in% rep_rebounds
  ) %>%
  mutate(
    `Mean delay (months)` = round(Mean_delay_months, 1),
    `Annual VL tests per 1000` = round(Total_VL_per_year, 0),
    `VL tests saved (%)` = round(Total_VL_saved * 100, 1)
  ) %>%
  dplyr::select(
    Scenario,
    `Annual rebound rate` = annual_rebound_rate,
    `AB specificity` = AB_Specificity,
    `AB rebound delay (years)` = AB_rebound_delay,
    `Mean delay (months)`,
    `Annual VL tests per 1000`,
    `VL tests saved (%)`
  ) %>%
  arrange(
    Scenario,
    `Annual rebound rate`,
    `AB specificity`,
    `AB rebound delay (years)`
  )

results_rep

main_text_table <- results %>%
  mutate(
    AB_Specificity = round(AB_Specificity, 2),
    AB_rebound_delay = round(AB_rebound_delay, 2),
    annual_rebound_rate = round(annual_rebound_rate, 2)
  ) %>%
  filter(
    Scenario %in% c("Annual AB", "biannual AB"),
    AB_Specificity %in% c(0.90, 0.95),
    AB_rebound_delay == 0.10,
    annual_rebound_rate %in% c(0.05, 0.10, 0.20)
  ) %>%
  mutate(
    `Mean delay (months)` = round(Mean_delay_months, 1),
    `Annual VL tests per 1000` = round(Total_VL_per_year, 0),
    `VL tests saved (%)` = round(Total_VL_saved * 100, 1)
  ) %>%
  select(
    Scenario,
    `Annual rebound rate` = annual_rebound_rate,
    `AB specificity` = AB_Specificity,
    `AB rebound delay (years)` = AB_rebound_delay,
    `Mean delay (months)`,
    `Annual VL tests per 1000`,
    `VL tests saved (%)`
  ) %>%
  arrange(Scenario, `Annual rebound rate`, `AB specificity`)

write.csv(main_text_table, 'table2_article.csv')

# knitr::kable(
#   main_text_table,
#   caption = "Representative main-text scenarios for antibody-guided monitoring."
# )

# =====================================================================
# UNIFIED PUBLICATION TYPOGRAPHY SETTINGS 
# =====================================================================
font_family     <- "sans"
axis_font_size  <- 24  # Clean, readable labels for 2D plots
tick_font_size  <- 20  # Visible numbers for 2D axis ticks

# --- BALANCED 3D FONTS ---
font_size_3d_axis <- 24  
font_size_3d_tick <- 18  

# =====================================================================
# PART 1: COMPILING & SAVING INDIVIDUAL 3D PLOTS (VIA WEBSHOT2)
# =====================================================================

# --- 1. Annual AB 3D Plot (p1_3d) ---
x_annual <- results %>% filter(Scenario == 'Annual AB')

p1_3d <- plot_ly(data = x_annual, x = ~AB_Specificity, y = ~annual_rebound_rate, z = ~Total_VL_saved,
                 type = "scatter3d", mode = "markers", 
                 marker = list(size = 5, color = ~Total_VL_saved, colorscale = "Viridis")) %>%
  layout(
    scene = list(
      xaxis = list(title = list(text = "AB Specificity", font = list(size = font_size_3d_axis, family = font_family)), tickfont = list(size = font_size_3d_tick)),
      yaxis = list(title = list(text = "Annual Rebound", font = list(size = font_size_3d_axis, family = font_family)), tickfont = list(size = font_size_3d_tick)),
      zaxis = list(title = list(text = "Total VL Saved", font = list(size = font_size_3d_axis, family = font_family)), tickfont = list(size = font_size_3d_tick),
                   range = c(0, max(x_annual$Total_VL_saved, na.rm = TRUE))),
      camera = list(eye = list(x = 1.6, y = 1.6, z = 1.1)), 
      aspectratio = list(x = 1.2, y = 1.2, z = 0.9) 
    ),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )

# Render and Save Individual Annual 3D Image
htmlwidgets::saveWidget(p1_3d, "temp_p1.html", selfcontained = TRUE)
webshot2::webshot("temp_p1.html", "Plot_1_Annual_3D.png", vwidth = 1200, vheight = 1000, zoom = 3.0)


# --- 2. Biannual AB 3D Plot (p2_3d) ---
x_biannual <- results %>% filter(Scenario == 'biannual AB')

p2_3d <- plot_ly(data = x_biannual, x = ~AB_Specificity, y = ~annual_rebound_rate, z = ~Total_VL_saved,
                 type = "scatter3d", mode = "markers", 
                 marker = list(size = 5, color = ~Total_VL_saved, colorscale = "Viridis")) %>%
  layout(
    scene = list(
      xaxis = list(title = list(text = "AB Specificity", font = list(size = font_size_3d_axis, family = font_family)), tickfont = list(size = font_size_3d_tick)),
      yaxis = list(title = list(text = "Annual Rebound", font = list(size = font_size_3d_axis, family = font_family)), tickfont = list(size = font_size_3d_tick)),
      zaxis = list(title = list(text = "Total VL Saved", font = list(size = font_size_3d_axis, family = font_family)), tickfont = list(size = font_size_3d_tick),
                   range = c(0, max(x_biannual$Total_VL_saved, na.rm = TRUE))),
      camera = list(eye = list(x = 1.6, y = 1.6, z = 1.1)),
      aspectratio = list(x = 1.2, y = 1.2, z = 0.9)
    ),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )

# Render and Save Individual Biannual 3D Image
htmlwidgets::saveWidget(p2_3d, "temp_p2.html", selfcontained = TRUE)
webshot2::webshot("temp_p2.html", "Plot_2_Biannual_3D.png", vwidth = 1200, vheight = 1000, zoom = 3.0)


# =====================================================================
# PART 2: COMPILING & SAVING INDIVIDUAL 2D PLOTS (VIA GGSAVE)
# =====================================================================

# --- 3. Annual Delay 2D Plot (gg_p3) ---
gg_p3 <- ggplot(results %>% filter(Scenario == 'Annual AB' & AB_Specificity == 0.6 & annual_rebound_rate == 0.04), 
                aes(y = Mean_delay_months, x = AB_rebound_delay)) + 
  geom_point(alpha = 0.8, size = 6, color = "#2b5c8f") +
  labs(x = "AB Rebound Delay", y = "Mean Delay (Months)") +
  theme_minimal(base_size = axis_font_size, base_family = font_family) +
  theme(
    axis.line        = element_line(colour = "black", linewidth = 0.8),
    axis.title.x     = element_text(face = "bold", margin = margin(t = 15)),
    axis.title.y     = element_text(face = "bold", margin = margin(r = 15)),
    axis.text        = element_text(size = tick_font_size, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92")
  )

# Save Individual Annual 2D Image
ggsave("Plot_3_Annual_2D.png", plot = gg_p3, width = 10, height = 8.33, dpi = 300)


# --- 4. Biannual Delay 2D Plot (gg_p4) ---
gg_p4 <- ggplot(results %>% filter(Scenario == 'biannual AB' & AB_Specificity == 0.6 & annual_rebound_rate == 0.04), 
                aes(y = Mean_delay_months, x = AB_rebound_delay)) + 
  geom_point(alpha = 0.8, size = 6, color = "#d95f02") +
  labs(x = "AB Rebound Delay", y = "Mean Delay (Months)") +
  theme_minimal(base_size = axis_font_size, base_family = font_family) +
  theme(
    axis.line        = element_line(colour = "black", linewidth = 0.8),
    axis.title.x     = element_text(face = "bold", margin = margin(t = 15)),
    axis.title.y     = element_text(face = "bold", margin = margin(r = 15)),
    axis.text        = element_text(size = tick_font_size, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92")
  )

# Save Individual Biannual 2D Image
ggsave("Plot_4_Biannual_2D.png", plot = gg_p4, width = 10, height = 8.33, dpi = 300)


# =====================================================================
# FILE SYSTEM CLEANUP & NOTIFICATION
# =====================================================================
if (file.exists("temp_p1.html")) file.remove("temp_p1.html")
if (file.exists("temp_p2.html")) file.remove("temp_p2.html")

cat("\nAll 4 files saved individually and completely independent of titles/composites!\n")


wide_rep_table <- results %>%
  mutate(
    AB_Specificity = round(AB_Specificity, 2),
    AB_rebound_delay = round(AB_rebound_delay, 2),
    annual_rebound_rate = round(annual_rebound_rate, 2),
    delay_label = paste0("Delay_", AB_rebound_delay)
  ) %>%
  filter(
    Scenario %in% c("Annual AB", "Biannual AB"),
    AB_Specificity %in% c(0.80, 0.90, 0.95),
    AB_rebound_delay %in% c(0.00, 0.10, 0.20),
    annual_rebound_rate %in% c(0.05, 0.10, 0.20)
  ) %>%
  mutate(
    metric_string = paste0(
      round(Mean_delay_months, 1), " mo; ",
      round(Total_VL_per_year, 0), " VL; ",
      round(Total_VL_saved * 100, 1), "% saved"
    )
  ) %>%
  select(
    Scenario,
    annual_rebound_rate,
    AB_Specificity,
    delay_label,
    metric_string
  ) %>%
  pivot_wider(
    names_from = delay_label,
    values_from = metric_string
  ) %>%
  arrange(Scenario, annual_rebound_rate, AB_Specificity)

knitr::kable(
  wide_rep_table,
  caption = "Representative scenarios in wide format. Entries are mean delay (months), annual VL tests per 1000, and VL tests saved (%)."
)

table_S1 <- results %>%
  filter(
    Scenario %in% c("Annual AB", "biannual AB"),
    AB_Specificity %in% c(0.80, 0.90, 0.95),
    AB_rebound_delay %in% c(0.0, 0.1, 0.2),
    annual_rebound_rate %in% c(0.05, 0.10, 0.20)
  ) %>%
  mutate(
    Mean_delay_months = round(Mean_delay_months, 1),
    VL_tests_per_1000 = round(Total_VL_per_year, 0),
    VL_saved_percent = round(Total_VL_saved * 100, 1)
  ) %>%
  select(
    Scenario,
    annual_rebound_rate,
    AB_Specificity,
    AB_rebound_delay,
    Mean_delay_months,
    VL_tests_per_1000,
    VL_saved_percent
  ) %>%
  arrange(Scenario, annual_rebound_rate, AB_Specificity)

write.csv(table_S1, "Table_S1_full_results.csv", row.names = FALSE)

table2_values <- results %>%
  filter(
    Scenario %in% c("Annual AB", "biannual AB"),
    AB_rebound_delay == 0.1,
    AB_Specificity %in% c(0.90, 0.95),
    annual_rebound_rate %in% c(0.10, 0.20)
  )
plot_data <- results %>%
  filter(
    Scenario == "Annual AB",
    AB_rebound_delay == 0.1
  ) %>%
  distinct(AB_Specificity, annual_rebound_rate, .keep_all = TRUE)

contour_plot <- ggplot(plot_data,
       aes(x = AB_Specificity,
           y = annual_rebound_rate,
           z = Total_VL_saved)) +
  
  geom_tile(aes(fill = Total_VL_saved)) +
  
  geom_contour(color = "white") +
  
  geom_text(
    stat = "contour",
    aes(label = after_stat(level)),
    size = 4
  ) +
  
  scale_fill_viridis_c(labels = scales::percent_format(accuracy = 1)) +
  
  labs(
    x = "AB Specificity",
    y = "Annual Rebound Rate",
    fill = "VL tests saved (%)"#,
    # title = "Annual Antibody Monitoring: VL Test Savings"
  ) +
  geom_point(
    data = table2_values,
    aes(
      x = ifelse(Scenario == "Annual AB", AB_Specificity - 0.005, AB_Specificity + 0.005),
      y = annual_rebound_rate,
      shape = Scenario,
      color = Scenario
    ),
    size = 3
  )
ggsave("contour_plot.png", plot = contour_plot, width = 13, height = 10, dpi = 300)
