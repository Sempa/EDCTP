library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(plotly)

# --- Parameter grid ---
spec_values <- seq(0.60, 0.95, 0.1)
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
