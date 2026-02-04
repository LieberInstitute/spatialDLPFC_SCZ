# Load libraries ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggplot2)
  library(rrvgo)
  library(sessioninfo)
})


# Load data ----
raw_gsea_prs_all <- read.csv(
  here(
    "processed-data/rds/14_prs_deg",
    "gsea_norm_PRS_DEG.csv"
  )
)

## Prep to format for rrvgo ----
gsea_prs_all <- raw_gsea_prs_all |> with(
  setNames(-log10(qvalue), ID)
)

# Prepare for rrvgo ----
## save semdata object ----
# NOTE: run one time to create the semdata object
# obj_semdata <- GOSemSim::godata(
#   annoDb = "org.Hs.eg.db",
#   ont = "BP" # ,
#   # keytype = keytype
# )

# saveRDS(
#   obj_semdata,
#   file = here(
#     "code/analysis/10_dx_deg_adjust_spd/01_GO_analysis",
#     "obj_semdata.rds"
#   )
# )

obj_semdata <- readRDS(
  here(
    "code/analysis/10_dx_deg_adjust_spd/01_GO_analysis",
    "obj_semdata.rds"
  )
)


# All terms together ----
## Calcualte and reduce the similarity matrix ----
simMatrix <- calculateSimMatrix(
  # The GO ID
  x = names(gsea_prs_all),
  orgdb = "org.Hs.eg.db",
  ont = "BP",
  method = "Rel",
  semdata = obj_semdata,
)

reducedTerms <- reduceSimMatrix(
  simMatrix,
  gsea_prs_all,
  threshold = 0.7,
  orgdb = "org.Hs.eg.db"
)

# Create treemap ----
pdf(
  file = here(
    "plots/14_prs_deg",
    "treemap_GO_prs_deg_layer_adj.pdf"
  )
)
treemapPlot(reducedTerms)
dev.off()
