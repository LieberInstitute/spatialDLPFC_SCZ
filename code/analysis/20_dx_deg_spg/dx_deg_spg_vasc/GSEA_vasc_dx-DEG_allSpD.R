# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(sessioninfo)
})

# Load Data ----
## Claudin-5+ SCZ-DEG results (all SpDs) ----
dx_res <- read_csv(
  here(
    "processed-data/20_dx_deg_spg/dx_deg_spg_vasc",
    "vasc-dx_DEG-allSpD.csv"
  )
)

# GSEA ----
## clusterProfiler ----
# Rank all genes by the moderated t-statistic of SCZ vs NTC
geneList <- dx_res |>
  arrange(desc(t_stat_scz)) |>
  dplyr::select(ensembl, t_stat_scz) |>
  deframe()

set.seed(20261007)
ego_gsea <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "BP", # focus on BP
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE,
  keyType = "ENSEMBL",
  pAdjustMethod = "fdr"
)

nrow(ego_gsea@result)

## Number of enriched terms by direction ----
ego_gsea@result |>
  count(direction = if_else(NES > 0, "up", "down"))

## Convert ENSMBLE IDs to SYMBOLs ----
ret_ego <- ego_gsea
ret_ego@result <- ret_ego@result |>
  mutate(
    core_enrichment_symbol = core_enrichment |> sapply(
      FUN = function(x) {
        str_split(x, "/")[[1]] |>
          bitr(fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db") |>
          pull(SYMBOL) |>
          paste0(collapse = "/")
      }
    )
  )

# Save results ----
write_csv(
  ret_ego@result,
  here(
    "processed-data/20_dx_deg_spg/dx_deg_spg_vasc",
    "gsea_vasc-dx_DEG-allSpD.csv"
  )
)

saveRDS(
  ret_ego,
  here(
    "processed-data/20_dx_deg_spg/dx_deg_spg_vasc",
    "gsea_vasc-dx_DEG-allSpD.rds"
  )
)

# Session Info ----
session_info()
