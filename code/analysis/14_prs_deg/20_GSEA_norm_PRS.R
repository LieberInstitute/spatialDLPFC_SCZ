# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  # library(DOSE)
  library(msigdbr)
  # library(fgsea)
  library(sessioninfo)
})

# Load Data ----
## PRS-deg gene list ----
test_res <- read_csv(
  file = here(
    "processed-data/rds/14_prs_deg/norm_PRS_deg",
    "norm_PRS_DEG_test_res_PRECAST07_donor_spd.csv"
  )
)

## Intro to MSigDB ----
# See all collection in MSigbr
# https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
msigdbr_collections()

# GSEA ----
## clusterProfiler ----
geneList <- test_res |>
  arrange(desc(t)) |>
  select(gene_id, t) |>
  deframe()

set.seed(20250723)
ego_gsea <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "BP", # focus on BP
  # ont = "ALL", # to have all
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE,
  keyType = "ENSEMBL",
  pAdjustMethod = "fdr"
)

nrow(ego_gsea@result)
# [1] 164

## Convert ENSMBLE IDs to SYMBOLs ----
ret_ego <- ego_gsea
ret_ego@result <- ret_ego@result |>
  mutate(
    core_enrichment_symbol = core_enrichment |> sapply(
      FUN = function(x) {
        # browser()
        str_split(x, "/")[[1]] |>
          bitr(fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db") |>
          pull(SYMBOL) |>
          paste0(collapse = "/")
      }
    )
  )

write_csv(
  ret_ego@result,
  here(
    "processed-data/rds/14_prs_deg",
    "gsea_norm_PRS_DEG.csv"
  )
)

# Session Info ----
session_info()
