# app.R
library(shiny)
library(gt)
library(dplyr)
library(gt)
library(scales)

ui <- fluidPage(
  titlePanel("ART Lab monitoring scenarios table"),
  sidebarLayout(
    sidebarPanel(
      numericInput("population", "Population (N):", value = 1000, min = 1, step = 1),
      sliderInput("AB_Specificity", "AB Specificity:", min = 0.5, max = 1, value = 0.95, step = 0.01),
      sliderInput("AB_rebound_delay", "AB rebound delay (years):", min = 0, max = 2, value = 0.2, step = 0.01),
      sliderInput("annual_rebound_rate", "Annual rebound rate (proportion):", min = 0, max = 1, value = 0.1, step = 0.01),
      width = 3
    ),
    mainPanel(
      # h4("Styled table (static, Excel-like layout)"),
      gt_output("table_gt"),
      br()#,
      # p("Note: calculations replicate the R code you provided. 'mean_delay_yrs' follows the same formula you supplied (AB_rebound_delay appears twice as in your snippet).")
    )
  )
)

server <- function(input, output, session) {
  
  tbl_calc <- reactive({
    
    population <- input$population
    AB_Specificity <- input$AB_Specificity
    AB_rebound_delay <- input$AB_rebound_delay
    annual_rebound_rate <- input$annual_rebound_rate
    
    monitoring_scenario <- c('Annual VL', 'Annual AB', 'Biannual AB')
    N <- rep(population, length(monitoring_scenario))
    prim_VL_tests <- c(population, 0, 0)
    testing_interval <- c(1, 1, 0.5)
    AB_testing <- c(0, rep(population,2))
    number_of_rebounds_PerAnnum <- rep(population * annual_rebound_rate, 3)
    number_of_rebounds_Per_test_round <- c(
      NA,
      population * annual_rebound_rate * (testing_interval[2] + AB_rebound_delay),
      (testing_interval[3] * annual_rebound_rate * population) + (population * annual_rebound_rate * AB_rebound_delay)
    )
    true_positives <- c(
      population[1]*annual_rebound_rate,
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
    false_negatives <- c(
      0,
      number_of_rebounds_Per_test_round[2] - true_positives[2],
      number_of_rebounds_Per_test_round[3] - true_positives[3]
    )
    true_negatives <- c(
      total_VL[1] - true_positives[1],
      (N[2] - number_of_rebounds_Per_test_round[2]) - false_positives[2],
      (N[3] - number_of_rebounds_Per_test_round[3]) - false_positives[3]
    )
    sensitivity <- c(
      1,
      testing_interval[2] / (testing_interval[2] + AB_rebound_delay),
      testing_interval[3] / (testing_interval[3] + AB_rebound_delay)
    )
    specificity <- c(
      1,
      AB_Specificity,
      AB_Specificity
    )
    positive_predictive_value <- c(
      1,
      # guard against division by zero
      ifelse((true_positives[2] + false_positives[2]) == 0, NA, true_positives[2] / (true_positives[2] + false_positives[2])),
      ifelse((true_positives[3] + false_positives[3]) == 0, NA, true_positives[3] / (true_positives[3] + false_positives[3]))
    )
    negative_predictive_value <- c(
      1,
      ifelse((true_negatives[2] + false_negatives[2]) == 0, NA, true_negatives[2] / (true_negatives[2] + false_negatives[2])),
      ifelse((true_negatives[3] + false_negatives[3]) == 0, NA, true_negatives[3] / (true_negatives[3] + false_negatives[3]))
    )
    total_VL_2 <- c(
      N[1],
      total_VL[2] / testing_interval[2],
      total_VL[3] / testing_interval[3]
    )
    total_AB <- c(
      0,
      AB_testing[2] / testing_interval[2],
      AB_testing[3] / testing_interval[3]
    )
    mean_delay_yrs <- c(
      0.5,
      0.5 * (AB_rebound_delay + testing_interval[2] + AB_rebound_delay), # kept exactly as supplied
      0.5 * (AB_rebound_delay + testing_interval[3] + AB_rebound_delay)
    )
    mean_delay_months <- mean_delay_yrs * 12
    
    # Assemble data.frame with columns that map to the Excel file column names.
    df <- tibble(
      "Scenario" = monitoring_scenario,
      "N" = N,
      "Prim_VL_tests" = prim_VL_tests,
      "Testing_interval_yrs" = testing_interval,
      "AB_testing_per_cycle" = AB_testing,
      "Rebounds_per_year" = number_of_rebounds_PerAnnum,
      "Rebounds_per_test_round" = number_of_rebounds_Per_test_round,
      "True_positives" = true_positives,
      "False_positives" = false_positives,
      "VL_confirmation" = vl_confirmation,
      "Total_VL" = total_VL,
      "False_negatives" = false_negatives,
      "True_negatives" = true_negatives,
      "Sensitivity" = sensitivity,
      "Specificity" = specificity,
      "PPV" = positive_predictive_value,
      "NPV" = negative_predictive_value,
      "Total_VL_per_yr" = total_VL_2,
      "Total_AB_per_yr" = total_AB,
      "Mean_delay_months" = mean_delay_months
    )
    
    # Format numeric columns nicely (rounding)
    df <- df %>%
      mutate(
        across(c(N, Prim_VL_tests, AB_testing_per_cycle, Rebounds_per_year, Rebounds_per_test_round,
                 True_positives, False_positives, VL_confirmation, Total_VL, False_negatives,
                 True_negatives, Total_VL_per_yr, Total_AB_per_yr),
               ~ round(.x, 1)),
        across(c(Sensitivity, Specificity, PPV, NPV), ~ round(.x, 3)),
        Mean_delay_months = round(Mean_delay_months, 1)
      )
    
    df
  })
  
  output$table_gt <- render_gt({
    df <- tbl_calc()
    
    # color palette for groups
    per_cycle_col <- "#E8F4FD"         # light blue
    per_test_round_col <- "#FFF4E6"    # light orange
    primary_perf_col <- "#E8FDE8"      # light green
    per_year_col <- "#F3E8FD"          # light purple
    
    # Build gt table
    gt_tbl <- df %>%
      gt(rowname_col = "Scenario") %>%
      cols_label(
        N = "N",
        Prim_VL_tests = "Prim VL tests",
        Testing_interval_yrs = "Testing interval (yrs)",
        AB_testing_per_cycle = "AB tests (per cycle)",
        Rebounds_per_year = "N rebounds P.A.",
        Rebounds_per_test_round = "N rebnd Per test round",
        True_positives = "True positives",
        False_positives = "False positives",
        VL_confirmation = "VL conf",
        Total_VL = "Tot VL",
        False_negatives = "False negatives",
        True_negatives = "True negatives",
        Sensitivity = "Sensitivity",
        Specificity = "Specificity",
        PPV = "PPV",
        NPV = "NPV",
        Total_VL_per_yr = "Tot VL (per yr)",
        Total_AB_per_yr = "Tot AB (per yr)",
        Mean_delay_months = "Mean delay (months)"
      ) %>%
      # Spanners (grouped headings)
      cols_spanner(
        label = "per cycle",
        columns = vars(N, Prim_VL_tests, Testing_interval_yrs, AB_testing_per_cycle)
      ) %>%
      cols_spanner(
        label = "per test round",
        columns = vars(Rebounds_per_test_round, True_positives, False_positives, VL_confirmation, Total_VL, False_negatives, True_negatives)
      ) %>%
      cols_spanner(
        label = "Primary test performance",
        columns = vars(Sensitivity, Specificity, PPV, NPV)
      ) %>%
      cols_spanner(
        label = "per year",
        columns = vars(Rebounds_per_year, Total_VL_per_yr, Total_AB_per_yr, Mean_delay_months)
      ) %>%
      # Styling group column backgrounds
      tab_style(
        style = list(cell_fill(color = per_cycle_col)),
        locations = cells_body(columns = c(N, Prim_VL_tests, Testing_interval_yrs, AB_testing_per_cycle))
      ) %>%
      tab_style(
        style = list(cell_fill(color = per_test_round_col)),
        locations = cells_body(columns = c(Rebounds_per_test_round, True_positives, False_positives, VL_confirmation, Total_VL, False_negatives, True_negatives))
      ) %>%
      tab_style(
        style = list(cell_fill(color = primary_perf_col)),
        locations = cells_body(columns = c(Sensitivity, Specificity, PPV, NPV))
      ) %>%
      tab_style(
        style = list(cell_fill(color = per_year_col)),
        locations = cells_body(columns = c(Rebounds_per_year, Total_VL_per_yr, Total_AB_per_yr, Mean_delay_months))
      ) %>%
      # Header formatting and numbers
      tab_header(
        title = md("Monitoring scenarios — computed table"),
        subtitle = md("Parameters: AB Specificity, AB rebound delay, Annual rebound rate")
      ) %>%
      fmt_number(columns = c(N, Prim_VL_tests, AB_testing_per_cycle, Rebounds_per_year, Rebounds_per_test_round,
                             True_positives, False_positives, VL_confirmation, Total_VL, False_negatives,
                             True_negatives, Total_VL_per_yr, Total_AB_per_yr),
                 decimals = 1) %>%
      fmt_number(columns = c(Sensitivity, Specificity, PPV, NPV), decimals = 3) %>%
      fmt_number(columns = c(Mean_delay_months), decimals = 1) %>%
      tab_options(
        table.font.size = px(14),
        row_group.font.size = px(14)
      )
    
    gt_tbl
  })
}

shinyApp(ui, server)
