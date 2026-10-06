# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(sessioninfo)
})

# Load Data ----
## Neuropil PRS-DEG results (GM only) ----
test_res <- read_csv(
  file = here(
    "processed-data/rds/14_prs_deg/Neuropil_PRS_DEG",
    "Neuropil_PRS_DEG_res_PRECAST07_donor_spd.csv"
  )
)

# GSEA ----
## clusterProfiler ----
# Rank all genes by the moderated t-statistic of norm_PRS
geneList <- test_res |>
  arrange(desc(t)) |>
  select(gene_id, t) |>
  deframe()

set.seed(20250723)
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
# [1] 444

## Number of enriched terms by direction ----
ego_gsea@result |>
  count(direction = if_else(NES > 0, "up", "down"))
#   direction   n
# 1      down 165
# 2        up 279

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
    "processed-data/rds/14_prs_deg/Neuropil_PRS_DEG",
    "gsea_Neuropil_PRS_DEG_GM.csv"
  )
)

saveRDS(
  ret_ego,
  here(
    "processed-data/rds/14_prs_deg/Neuropil_PRS_DEG",
    "gsea_Neuropil_PRS_DEG_GM.rds"
  )
)

# Session Info ----
session_info()
