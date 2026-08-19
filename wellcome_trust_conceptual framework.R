install.packages("DiagrammeR")
install.packages("DiagrammeRsvg")
install.packages("rsvg")

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

g <- grViz("
digraph conceptual_framework {

  graph [layout = dot, rankdir = TB]

  node [shape = rectangle, style = filled, fillcolor = white, fontsize = 12]

  A [label = 'HIV Viral Dynamics\n(suppression ↔ breakthrough)']
  B [label = 'Antigen Exposure Profile\n(time-varying + cumulative viraemia)']
  C [label = 'Antibody Immune Response\n• Contraction (λ)\n• Expansion (r)\n• Avidity maturation']
  D [label = 'Measured Antibody Trajectories\n(venous | capillary | oral fluid)']

  M [label = 'Modifying Factors:\n• HIV subtype\n• Resistance vs non-resistance\n• Treatment & regimen\n• Host factors',
     fillcolor = lightgrey]

  A -> B -> C -> D

  { rank = same; C; M }

  M -> C [style = dashed]
}
")

# Convert diagram to SVG
svg <- export_svg(g)

# Write to PNG
rsvg_png(charToRaw(svg), file = "conceptual_framework.png")


##################################################################################
##Gantt Chart
#################################################################################

library(ggplot2)

# Create dataset (time in years; 0–0.25 = first 3 months)
gantt <- data.frame(
  
  Task = c(
    "Ethics approval",
    
    "Cohort harmonisation",
    "Assay standardisation",
    "Initial kinetic modelling",
    
    "Subtype comparative analysis",
    "Resistance pathway modelling",
    
    "Compartment data collection",
    "Compartment modelling",
    
    "Integrated host–virus model"
  ),
  
  Start = c(
    0,
    0.25, 0.25, 0.25,
    2, 2,
    3, 3,
    4
  ),
  
  End = c(
    0.25,
    2, 2, 2,
    3, 3,
    5, 5,
    5
  ),
  
  Aim = c(
    "Setup",
    
    "Aim 1", "Aim 1", "Aim 1",
    
    "Aim 1", "Aim 2",
    
    "Aim 3", "Aim 3",
    
    "Integration"
  )
)

# Order tasks neatly
gantt$Task <- factor(gantt$Task, levels = rev(gantt$Task))

# Plot
p <- ggplot(gantt, aes(x = Start, xend = End, y = Task, yend = Task, color = Aim)) +
  
  geom_segment(size = 6) +
  
  # Milestones
  geom_vline(xintercept = c(0.25, 2, 3, 5), linetype = "dashed") +
  
  scale_x_continuous(
    breaks = c(0, 0.25, 1, 2, 3, 4, 5),
    labels = c("Start", "3 mo", "Year 1", "Year 2", "Year 3", "Year 4", "Year 5")
  ) +
  
  scale_color_manual(values = c(
    "Setup" = "black",
    "Aim 1" = "steelblue",
    "Aim 2" = "darkorange",
    "Aim 3" = "forestgreen",
    "Integration" = "purple"
  )) +
  
  labs(
    title = "Research Plan Timeline and Milestones",
    subtitle = "Staged delivery aligned to specific aims",
    x = "Project Timeline",
    y = "",
    color = "Project Phase"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )

print(p)

# Save high-resolution image
ggsave("gantt_chart_wellcome.png", width = 11, height = 6, dpi = 300)
