# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggplot2)
  library(readxl)
  library(rrvgo)
  library(sessioninfo)
})

# Load data ----
## Annotated GO terms from Sang Ho ----
raw_go_gsea_172 <- read_xlsx(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd/SHK_annotation", "gsea_PRECAST07_SHK.xlsx"
  )
)

sum(raw_go_gsea_172$p.adjust < 0.05)
# NOTE: the file contains only stat sig GO terms at 0.05 level.

## Organize format for rrvgo ----
go_gsea_172 <- raw_go_gsea_172 |> with(
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
  x = names(go_gsea_172),
  orgdb = "org.Hs.eg.db",
  ont = "BP",
  method = "Rel",
  semdata = obj_semdata,
)

reducedTerms <- reduceSimMatrix(
  simMatrix,
  go_gsea_172,
  threshold = 0.7,
  orgdb = "org.Hs.eg.db"
)

# Create treemap ----
pdf(
  file = here(
    "plots/10_dx_deg_adjust_spd",
    "treemap_GO_layer_restrict.pdf"
  )
)
treemapPlot(reducedTerms)
dev.off()


## up-regulated terms (NES > 0) ----
go_gsea_172_up <- raw_go_gsea_172 |>
  filter(NES > 0) |>
  with(
    setNames(-log10(qvalue), ID)
  )

simMatrix_up <- calculateSimMatrix(
  # The GO ID
  x = names(go_gsea_172_up),
  orgdb = "org.Hs.eg.db",
  ont = "BP",
  method = "Rel",
  semdata = obj_semdata,
)

reducedTerms_up <- reduceSimMatrix(
  simMatrix_up,
  go_gsea_172_up,
  threshold = 0.7,
  orgdb = "org.Hs.eg.db"
)

# Create treemap ----
pdf(
  file = here(
    "plots/10_dx_deg_adjust_spd/archived",
    "treemap_GO_layer_restrict.pdf"
  )
)
treemapPlot(reducedTerms_up)
dev.off()

## Down-reg terms (NES < 0) ----
go_gsea_172_down <- raw_go_gsea_172 |>
  filter(NES < 0) |>
  with(
    setNames(-log10(qvalue), ID)
  )

simMatrix_down <- calculateSimMatrix(
  # The GO ID
  x = names(go_gsea_172_down),
  orgdb = "org.Hs.eg.db",
  ont = "BP",
  method = "Rel",
  semdata = obj_semdata,
)

reducedTerms_down <- reduceSimMatrix(
  simMatrix_down,
  go_gsea_172_down,
  threshold = 0.7,
  orgdb = "org.Hs.eg.db"
)

pdf(
  file = here(
    "plots/10_dx_deg_adjust_spd/archived",
    "treemap_GO_layer_restrict.pdf"
  )
)
treemapPlot(reducedTerms_down)
dev.off()


# Session info ----
session_info()
