library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

graph_code <- "
digraph monitoring_pathways {

graph [
layout = dot,
rankdir = LR
]

node [
shape = box,
style = 'rounded,filled',
fillcolor = white,
fontname = Helvetica,
fontsize = 33,
penwidth = 1.8,
margin = 0.25
]

edge [
fontname = Helvetica,
fontsize = 20,
penwidth = 1.8
]

# -------------------------------------
# Annual VL pathway
# -------------------------------------

subgraph cluster_1 {

label = 'A. Annual VL Monitoring (Standard of Care)'
fontsize = 30
fontname = Helvetica
penwidth = 2.2

A1 [label='Patient on ART']
A2 [label='Annual VL Test']
A3 [label='Virological Rebound Detected']
A4 [label='Clinical Action']

A1 -> A2 -> A3 -> A4
}

# -------------------------------------
# Annual AB pathway
# -------------------------------------

subgraph cluster_2 {

label = 'B. Annual Antibody Monitoring'
fontsize = 30
fontname = Helvetica
penwidth = 2.2

B1 [label='Patient on ART']
B2 [label='Annual Antibody Test']
B3 [label='AB Positive?']

B4 [label='Continue ART']
B5 [label='Confirmatory VL Test']
B6 [label='Virological Rebound Detected']
B7 [label='Clinical Action']

B1 -> B2 -> B3
B3 -> B4 [label='No']
B3 -> B5 [label='Yes']
B5 -> B6 -> B7
}

# -------------------------------------
# Biannual AB pathway
# -------------------------------------

subgraph cluster_3 {

label = 'C. Biannual Antibody Monitoring'
fontsize = 30
fontname = Helvetica
penwidth = 2.2

C1 [label='Patient on ART']
C2 [label='Antibody Test Every 6 Months']
C3 [label='AB Positive?']

C4 [label='Continue ART']
C5 [label='Confirmatory VL Test']
C6 [label='Virological Rebound Detected']
C7 [label='Clinical Action']

C1 -> C2 -> C3
C3 -> C4 [label='No']
C3 -> C5 [label='Yes']
C5 -> C6 -> C7
}

}
"

# Create graph
g <- grViz(graph_code)

# Export as SVG
svg <- DiagrammeRsvg::export_svg(g)

# Ensure svg is a single character string
svg <- paste(svg, collapse = "\n")

# Save SVG
writeLines(svg, "Figure1_MonitoringPathways.svg")

# Export PNG
rsvg::rsvg_png(
  charToRaw(svg),
  file = "Figure1_MonitoringPathways.png",
  width = 4000,
  height = 2200
)
