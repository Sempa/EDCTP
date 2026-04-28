# cohort_model_schematic.R
# Generates a cohort-model schematic for:
#   Strategy A: PCR-only monitoring
#   Strategy B: Ab triage -> confirmatory PCR
#
# Output: cohort_model_schematic.svg and cohort_model_schematic.png
#
# Packages:
# install.packages(c("DiagrammeR","DiagrammeRsvg","rsvg"))

library(DiagrammeR)
library(rsvg)
library(glue)

# ---- Edit labels here if you want to customize text ----
labels <- list(
  
  cohort =
    "Cohort of adults receiving ART
(N persons entering simulation)",
  
  pcr_all =
    "Routine HIV RNA PCR monitoring
1 test per person-year
Cost = cPCR + cvisit",
  
  pcr_pos =
    "Detectable viral load
(VL ≥1000 copies/mL)",
  
  pcr_neg =
    "Viral suppression maintained",
  
  manage =
    "Clinical management after detection:
adherence support or regimen switch
based on cause of rebound",
  
  note_pcr =
    "Outputs:
Annual cost
Suppression-years
Life-years
DALYs",
  
  ab_triage =
    "Antibody triage at frequency f
(annual to 6-weekly)
Cost = cAb + cvisit",
  
  ab_pos =
    "Reactive antibody test
Proceed to confirmatory PCR",
  
  pcr_conf =
    "Confirmatory PCR
performed only after
positive antibody result",
  
  manage2 =
    "Earlier detection may improve:
re-suppression probability
reduce mortality
reduce disability",
  
  net =
    "Incremental outcomes:
Δ Cost
Δ Suppression-years
Δ Life-years
Δ DALYs

Primary ICER:
Cost per DALY averted"
)

# ---- Graphviz (DOT) specification ----
dot <- glue('
digraph cohort_model {{

graph [
  layout = dot,
  rankdir = TB,
  fontsize = 20,
  labelloc = "t",
  # label = "Figure 1. Decision-analytic cohort model comparing PCR-only monitoring with antibody triage",
  splines = true
]

node [
  shape = box,
  style = "rounded,filled",
  fontname = Helvetica,
  fontsize = 11,
  color = "#111827",
  penwidth = 1.2
]

edge [
  fontname = Helvetica,
  fontsize = 10,
  color = "#111827",
  penwidth = 1.1
]

subgraph cluster_pcr {{
  label = "A. PCR-only monitoring";
  color = "#c7d2fe";
  style = "rounded";

  cohort1   [label="{labels$cohort}", fillcolor="#eef2ff"];
  pcr_all   [label="{labels$pcr_all}", fillcolor="#ecfeff"];
  pcr_pos   [label="{labels$pcr_pos}", fillcolor="#fff7ed"];
  pcr_neg   [label="{labels$pcr_neg}", fillcolor="#f0fdf4"];
  manage1   [label="{labels$manage}", fillcolor="#fefce8"];
  outcomes1 [label="{labels$note_pcr}", fillcolor="#eef2ff"];

  cohort1 -> pcr_all;
  pcr_all -> pcr_pos [label="VL rebound"];
  pcr_all -> pcr_neg [label="No rebound"];
  pcr_pos -> manage1;
  manage1 -> outcomes1;
}}

subgraph cluster_triage {{
  label = "B. Antibody triage strategy";
  color = "#c7d2fe";
  style = "rounded";

  cohort2   [label="{labels$cohort}", fillcolor="#eef2ff"];
  ab_test   [label="{labels$ab_triage}", fillcolor="#ecfeff"];
  ab_pos    [label="{labels$ab_pos}", fillcolor="#fff7ed"];
  ab_neg    [label="No reactive Ab result", fillcolor="#f0fdf4"];
  pcr_conf  [label="{labels$pcr_conf}", fillcolor="#fefce8"];
  manage2   [label="{labels$manage2}", fillcolor="#ffffff"];
  outcomes2 [label="{labels$net}", fillcolor="#eef2ff"];

  cohort2 -> ab_test;
  ab_test -> ab_pos [label="Positive"];
  ab_test -> ab_neg [label="Negative"];
  ab_pos -> pcr_conf;
  pcr_conf -> manage2;
  manage2 -> outcomes2;
}}

{{rank = same; cohort1; cohort2}}

}}
')

# Render diagram in RStudio Viewer / notebook
viz <- grViz(dot)
print(viz)

# ---- Export to SVG and PNG (optional but recommended) ----
# These require DiagrammeRsvg + rsvg
if (requireNamespace("DiagrammeRsvg", quietly = TRUE)) {
  svg_txt <- DiagrammeRsvg::export_svg(viz)
  writeLines(svg_txt, "cohort_model_schematic.svg")
  
  if (requireNamespace("rsvg", quietly = TRUE)) {
    rsvg::rsvg_png(charToRaw(svg_txt), "cohort_model_schematic.png",
                   width = 2400, height = 1200)
  } else {
    message("Package 'rsvg' not installed: SVG saved, PNG not created.")
  }
} else {
  message("Package 'DiagrammeRsvg' not installed: diagram rendered but not exported.")
}

# Convert SVG file to PNG
rsvg_png("CE_AB_analysis/cohort_model_schematic.svg", "CE_AB_analysis/cohort_model_schematic.png",
         width = 2400, height = 1200)
