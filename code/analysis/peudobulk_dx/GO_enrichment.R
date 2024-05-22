# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  # library(DOSE)
  library(sessioninfo)
})

# Load Data ----
gene_df <- read_csv(
  #   here(
  #     "processed-data/PB_dx_genes/",
  #     "test_PRECAST_07.csv"
  #   )
  # )

  "~/Downloads/test_PRECAST_07.csv"
)

sig_gene_df <- gene_df |>
  filter(fdr_ntc <= 0.05)

sig_gene <- sig_gene_df |>
  pull(ensembl)


up_gene <- gene_df |>
  filter(
    fdr_ntc <= 0.05,
    logFC_scz > 0
  ) |>
  pull(ensembl)

down_gene <- gene_df |>
  filter(
    fdr_ntc <= 0.05,
    logFC_scz < 0
  ) |>
  pull(ensembl)

# Export the genes for website
# query_string <- paste0(sig_gene, collapse = "::")


dir.create(
  here("processed-data/PB_dx_genes/enrichment"),
  showWarnings = FALSE
)

# GO enrichment (Overrepresentative analysis) ----
gc_list <- list(
  "All" = sig_gene,
  "Up" = up_gene,
  "Down" = down_gene

)

GO_ora <- compareCluster(
  geneCluster = gc_list, fun = enrichGO,
  universe = gene_df$ensembl,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "ENSEMBL",
  pAdjustMethod = "BH",
  # pvalueCutoff = 0.01,
  # qvalueCutoff = 0.05,
  readable = TRUE
)


# Visualization ----
dotplot(GO_ora, showCategory = 10, split = "ONTOLOGY", includeAll = FALSE, group = TRUE)+
scale_x_discrete(limits = c("Up", "Down", "All"))

# Results ----
GO_ora@compareClusterResult |> View()





# Deprecated Code ----
## All genes ----
# ego <- enrichGO(
#   gene = sig_gene,
# universe = gene_df$ensembl,
# OrgDb = org.Hs.eg.db,
# ont = "ALL",
# keyType = "ENSEMBL",
# pAdjustMethod = "BH",
# pvalueCutoff = 0.01,
# qvalueCutoff = 0.05,
# readable = TRUE
# )

# # dotplot(ego)

# write.csv(
#   ego@result,
#   here(
#     "processed-data/PB_dx_genes/enrichment",
#     "GO_ERA_PRECAST_07.csv"
#   )
# )

# ## Up reg genes ----
# up_gene <- gene_df |>
#   filter(
#     fdr_ntc <= 0.05,
#     logFC_scz > 0
#   ) |>
#   pull(ensembl)
# length(up_gene)
# ego_up <- enrichGO(
#   gene = up_gene,
#   universe = gene_df$ensembl,
#   OrgDb = org.Hs.eg.db,
#   ont = "ALL",
#   keyType = "ENSEMBL",
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.01,
#   qvalueCutoff = 0.05,
#   readable = TRUE
# )
# write.csv(
#   ego_up@result,
#   here(
#     "processed-data/PB_dx_genes/enrichment",
#     "GO_ERA_up_gene_PRECAST_07.csv"
#   )
# )
# ## Down reg genes ----
# down_gene <- gene_df |>
#   filter(
#     fdr_ntc <= 0.05,
#     logFC_scz < 0
#   ) |>
#   pull(ensembl)
# length(down_gene)
# ego_down <- enrichGO(
#   gene = down_gene,
#   universe = gene_df$ensembl,
#   OrgDb = org.Hs.eg.db,
#   ont = "ALL",
#   keyType = "ENSEMBL",
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.01,
#   qvalueCutoff = 0.05,
#   readable = TRUE
# )
# write.csv(
#   ego_down@result,
#   here(
#     "processed-data/PB_dx_genes/enrichment",
#     "GO_ERA_down_gene_PRECAST_07.csv"
#   )
# )
# # Overlaping pathways between down and up regulated genes
# intersect(
#   ego_up@result$ID,
#   ego_down@result$ID
# )
# # 0
# diff_processes <- setdiff(
#   union(
#     ego_up@result$ID,
#     ego_down@result$ID
#   ),
#   ego@result$ID
# )
# # Which are the processes not in the lump sum analysis?
# rbind(
#     ego_up@result,
#     ego_down@result
#   ) |> filter(ID %in% diff_processes) |> View()

# Session Info ----
session_info()
