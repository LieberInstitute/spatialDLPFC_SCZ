# Load libraries ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(rrvgo)
  library(sessioninfo)
})

# Load data ----
raw_gsea_prs <- read_csv(
  here(
    "processed-data/rds/14_prs_deg/Neuropil_PRS_DEG",
    "gsea_Neuropil_PRS_DEG_GM.csv"
  )
)

## Semantic data for GO BP ----
# NOTE: created in code/analysis/14_prs_deg/21_treemap_GSEA_norm_PRS.R
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
  pdf(file = here("plots/14_prs_deg/Neuropil_PRS_DEG", file_name))
  treemapPlot(reducedTerms, title = title)
  dev.off()

  reducedTerms
}

# All terms together ----
reduced_all <- reduce_and_plot(
  raw_gsea_prs,
  "treemap_GO_Neuropil_PRS_DEG_GM.pdf",
  title = "Neuropil PRS-DEG GSEA (GO BP): all terms"
)

# By direction of enrichment ----
## Positively enriched (higher expression with higher PRS) ----
reduced_up <- reduce_and_plot(
  raw_gsea_prs |> filter(NES > 0),
  "treemap_GO_Neuropil_PRS_DEG_GM_up.pdf",
  title = "Neuropil PRS-DEG GSEA (GO BP): NES > 0"
)

## Negatively enriched (lower expression with higher PRS) ----
reduced_down <- reduce_and_plot(
  raw_gsea_prs |> filter(NES < 0),
  "treemap_GO_Neuropil_PRS_DEG_GM_down.pdf",
  title = "Neuropil PRS-DEG GSEA (GO BP): NES < 0"
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
      "processed-data/rds/14_prs_deg/Neuropil_PRS_DEG",
      "rrvgo_reduced_terms_Neuropil_PRS_DEG_GM.csv"
    )
  )

# Session info ----
session_info()
