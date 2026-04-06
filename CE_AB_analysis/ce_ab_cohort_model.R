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

# ---- Edit labels here if you want to customize text ----
labels <- list(
  cohort = "Cohort of PLWH monitored\n(N persons; annual probability of viraemia = p_v)",
  
  pcr_all = "Routine HIV RNA (PCR) for everyone\n(1 test/person-year; cost = c_PCR + c_visit)",
  pcr_pos = "PCR: VL ≥ 1000\n(True viraemia detected)",
  pcr_neg = "PCR: VL < 1000\n(Virally suppressed)",
  manage = "Clinical response (EAC ± regimen change)\nProbability of re-suppression: p_res\n(weighted by resistance vs behaviour)",
  resup = "Re-suppressed earlier\n(benefit: Δt years gained)",
  noresup = "Not re-suppressed\n(no benefit)",
  note_pcr = "Outcomes used in code: cost_pcronly, eff_pcronly\n(assumes average detection delay ≈ interval/2)",
  
  ab_triage = "Antibody (Ab) triage at chosen frequency\n(f tests/year; cost = c_Ab + c_visit)",
  ab_pos = "Ab positive\n(True pos: sens among viraemic\nFalse pos: 1−spec among suppressed)",
  ab_neg = "Ab negative\n(False neg possible if viraemic)",
  pcr_conf = "Confirmatory PCR for Ab+ only\n(expected PCR/person-year = f·P(Ab+)\n(additional cost = c_PCR + c_visit)",
  manage2 = "If PCR confirms VL ≥ 1000 → clinical response\nRe-suppression probability: p_res (weighted)\nEarlier detection benefit depends on Δdelay vs PCR-only",
  net = "Net outcomes used in code:\nΔcost = cost_triage − cost_pcronly;\nΔeff = eff_triage − eff_pcronly; ICER = Δcost/Δeff",
  note_triage = "Note: Ab− branch implies missed/late detection for some viraemia\n(in code simplified via sens and Δdelay)."
)

# ---- Graphviz (DOT) specification ----
dot <- sprintf('
  digraph cohort_model {
    graph [layout = dot, rankdir = TB, fontsize = 18, labelloc = "t",
           label = "",#Cohort model schematic corresponding to ce_ab_triage.R
           splines = true]
    node  [shape = box, style = "rounded,filled", fontname = Helvetica, color = "#111827",
           fontsize = 11, penwidth = 1.2]
    edge  [fontname = Helvetica, fontsize = 10, color = "#111827", penwidth = 1.1]

    // --- Left panel: PCR-only ---
    subgraph cluster_pcr {
      label = "Strategy A: PCR-only monitoring";
      color = "#c7d2fe";
      style = "rounded";

      cohort1 [label = "%s", fillcolor = "#eef2ff"];
      pcr_all [label = "%s", fillcolor = "#ecfeff"];
      pcr_pos [label = "%s", fillcolor = "#fff7ed"];
      pcr_neg [label = "%s", fillcolor = "#f0fdf4"];
      manage1 [label = "%s", fillcolor = "#fefce8"];
      resup   [label = "%s", fillcolor = "#f0fdf4"];
      noresup [label = "%s", fillcolor = "#fee2e2"];
      note1   [label = "%s", fillcolor = "#ffffff", fontsize = 9];

      cohort1 -> pcr_all;
      pcr_all -> pcr_pos [label = "p_v"];
      pcr_all -> pcr_neg [label = "1−p_v"];
      pcr_pos -> manage1;
      manage1 -> resup   [label = "p_res"];
      manage1 -> noresup [label = "1−p_res"];
      resup   -> note1   [style = invis];
      noresup -> note1   [style = invis];

      {rank = same; pcr_pos; pcr_neg}
      {rank = same; resup; noresup}
    }

    // --- Right panel: Ab triage -> PCR ---
    subgraph cluster_triage {
      label = "Strategy B: Ab triage → confirmatory PCR";
      color = "#c7d2fe";
      style = "rounded";

      cohort2 [label = "%s", fillcolor = "#eef2ff"];
      ab_tri  [label = "%s", fillcolor = "#ecfeff"];
      ab_pos  [label = "%s", fillcolor = "#fff7ed"];
      ab_neg  [label = "%s", fillcolor = "#f0fdf4"];
      pcr_conf [label = "%s", fillcolor = "#fefce8"];
      manage2  [label = "%s", fillcolor = "#ffffff"];
      net      [label = "%s", fillcolor = "#eef2ff"];
      note2    [label = "%s", fillcolor = "#ffffff", fontsize = 9];

      cohort2 -> ab_tri;
      ab_tri -> ab_pos [label = "P(Ab+) = p_v·sens + (1−p_v)·(1−spec)"];
      ab_tri -> ab_neg [label = "P(Ab−)"];
      ab_pos -> pcr_conf;
      pcr_conf -> manage2;
      manage2 -> net;
      net -> note2;

      {rank = same; ab_pos; ab_neg}
    }

    // Make both cohort nodes align at top
    {rank = same; cohort1; cohort2}
  }
',
               labels$cohort,
               labels$pcr_all,
               labels$pcr_pos,
               labels$pcr_neg,
               labels$manage,
               labels$resup,
               labels$noresup,
               labels$note_pcr,
               labels$cohort,
               labels$ab_triage,
               labels$ab_pos,
               labels$ab_neg,
               labels$pcr_conf,
               labels$manage2,
               labels$net,
               labels$note_triage
)

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
rsvg_png("cohort_model_schematic.svg", "cohort_model_schematic.png",
         width = 2400, height = 1200)