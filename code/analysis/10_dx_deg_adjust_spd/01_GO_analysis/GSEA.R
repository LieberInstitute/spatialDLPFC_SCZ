# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  # library(org.Hs.eg.db)
  # library(clusterProfiler)
  # library(DOSE)
  library(msigdbr)
  library(fgsea)
  library(sessioninfo)
})

# Load Data ----
## dx-deg gene list ----
gene_df <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_PRECAST07.csv"
  )
)


## Intro to MSigDB ----
# See all collection in MSigbr
# https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
msigdbr_collections()

# GSEA ----

## clusterProfiler ----
geneList <- gene_df |>
  arrange(desc(t_stat_scz)) |>
  select(ensembl, t_stat_scz) |>
  deframe()

set.seed(20250422)
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
# [1] 294


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
    "processed-data/rds/10_dx_deg_adjust_spd",
    "gsea_PRECAST07.csv"
  )
)


## fgsea (hard to format to human readable format) ----
### set seed for reproducibility ----
set.seed(20250421)

msigdbr_gene_sets <- msigdbr(species = "Homo sapiens", category = "C5")
msigdbr_list <- split(
  x = msigdbr_gene_sets$ensembl_gene,
  f = msigdbr_gene_sets$gs_name
)

up_rank_gene_list <- gene_df |>
  arrange(desc(t_stat_scz)) |>
  select(ensembl, t_stat_scz) |>
  deframe()

set.seed(20250421)
up_reg_gsea <- fgsea(
  pathways = msigdbr_list,
  stats = up_rank_gene_list,
  nPermSimple = 10000
)

up_reg_gsea |>
  arrange(padj) |>
  filter(padj <= 0.05) |>
  View()


# SYnapse_pruning is down regulated, concordance with the GO ORA results


tmp_collap <- collapsePathways(
  up_reg_gsea |> filter(padj <= 0.05),
  msigdbr_list,
  up_rank_gene_list
)


# Visualizaiton ---
plotGseaTable(
  msigdbr_list["GOBP_SYNAPSE_PRUNING"],
  up_rank_gene_list, up_reg_gsea,
  gseaParam = 0.5
)

plotEnrichment(
  msigdbr_list["GOBP_SYNAPSE_PRUNING"],
  up_rank_gene_list
) + labs(title = "Programmed Cell Death")
