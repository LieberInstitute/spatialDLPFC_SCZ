# Load library  ----
suppressPackageStartupMessages()({
  library(here)
  library(spatialLIBD)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  # library(limma)
  # library(readxl)
  library(sessioninfo)
})


# Load data ----
dx_res <- read_csv(
  here(
    "code/analysis/dx_deg_spg_neuropil",
    "neuropil-dx_DEG-GM.csv"
  )
)


# GO Enrichment ----
## Overrepresentation test ----
up_gene <- dx_res |>
  filter(
    logFC_scz > 0,
    fdr_scz < 0.10
  ) |>
  pull(ensembl)

down_gene <- dx_res |>
  filter(
    logFC_scz < 0,
    fdr_scz < 0.10
  ) |>
  pull(ensembl)

# Export the genes for website
# query_string <- paste0(sig_gene, collapse = "::")


# Up-reg genes -----
up_ego <- enrichGO(
  gene = up_gene,
  universe = dx_res$ensembl,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "ENSEMBL",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

up_ego@result |>
  write.csv(
    here(
      # TODO:
      # "processed-data/rds/10_dx_deg_adjust_spd",
      "code/analysis/dx_deg_spg_neuropil",
      "GO_neuropil_GW_up_gene.csv"
    )
  )

# Down-reg genes ----
down_ego <- enrichGO(
  gene = down_gene,
  universe = dx_res$ensembl,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "ENSEMBL",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

down_ego@result |>
  write.csv(
    here(
      # TODO:
      # "processed-data/rds/10_dx_deg_adjust_spd",
      "code/analysis/dx_deg_spg_neuropil",
      "GO_neuropil_GW_down_gene.csv"
    )
  )

# TODO:

## GSEA test ----
set.seed(20250428)

geneList <- dx_res |>
  arrange(desc(t_stat_scz)) |>
  dplyr::select(ensembl, t_stat_scz) |>
  deframe()

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
# [1] 355


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
      # TODO:
      # "processed-data/rds/10_dx_deg_adjust_spd",
      "code/analysis/dx_deg_spg_neuropil",
      "GSEA_neuropil_GW_down_gene.csv"
  )
)
