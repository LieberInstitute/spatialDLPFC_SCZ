# Load libraries ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(rrvgo)
  library(sessioninfo)
})

# Load data ----
# NOTE: created in code/analysis/20_dx_deg_spg/dx_deg_spg_pnn/GSEA_pnn_dx-DEG_GM.R
raw_gsea_dx <- read_csv(
  here(
    "processed-data/20_dx_deg_spg/dx_deg_spg_pnn",
    "gsea_pnn-dx_DEG-GM.csv"
  )
)

## Semantic data for GO BP ----
obj_semdata <- readRDS(
  here(
    "code/analysis/10_dx_deg_adjust_spd/01_GO_analysis",
    "obj_semdata.rds"
  )
)

# Helper: reduce GO terms and draw treemap ----
# gsea_df: GSEA result table with columns ID and qvalue
# returns the reducedTerms data.frame (NULL if too few terms)
reduce_and_plot <- function(gsea_df, file_name, title = "") {
  # Scores for rrvgo: higher is more important
  scores <- gsea_df |> with(setNames(-log10(qvalue), ID))

  if (length(scores) < 2) {
    message("Fewer than 2 GO terms for ", file_name, "; skip treemap.")
    return(NULL)
  }

  ## Calculate and reduce the similarity matrix ----
  simMatrix <- calculateSimMatrix(
    x = names(scores),
    orgdb = "org.Hs.eg.db",
    ont = "BP",
    method = "Rel",
    semdata = obj_semdata
  )

  reducedTerms <- reduceSimMatrix(
    simMatrix,
    scores,
    threshold = 0.7,
    orgdb = "org.Hs.eg.db"
  )

  ## Create treemap ----
  pdf(file = here("plots/20_dx_deg_spg/dx_deg_spg_pnn", file_name))
  treemapPlot(reducedTerms, title = title)
  dev.off()

  reducedTerms
}

# All terms together ----
reduced_all <- reduce_and_plot(
  raw_gsea_dx,
  "treemap_GO_pnn-dx_DEG-GM.pdf",
  title = "PNN+ SCZ-DEG GSEA (GO BP): all terms"
)

# By direction of enrichment ----
## Positively enriched (higher expression in SCZ) ----
reduced_up <- reduce_and_plot(
  raw_gsea_dx |> filter(NES > 0),
  "treemap_GO_pnn-dx_DEG-GM_up.pdf",
  title = "PNN+ SCZ-DEG GSEA (GO BP): NES > 0"
)

## Negatively enriched (lower expression in SCZ) ----
reduced_down <- reduce_and_plot(
  raw_gsea_dx |> filter(NES < 0),
  "treemap_GO_pnn-dx_DEG-GM_down.pdf",
  title = "PNN+ SCZ-DEG GSEA (GO BP): NES < 0"
)

# Save reduced terms (parent clusters) ----
bind_rows(
  all = reduced_all,
  up = reduced_up,
  down = reduced_down,
  .id = "set"
) |>
  write_csv(
    here(
      "processed-data/20_dx_deg_spg/dx_deg_spg_pnn",
      "rrvgo_reduced_terms_pnn-dx_DEG-GM.csv"
    )
  )

# Session info ----
session_info()
