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
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_PRECAST07.csv"
  )
)
adj_p_cutoff <- 0.10

sig_gene_df <- gene_df |>
  filter(fdr_ntc <= adj_p_cutoff)

sig_gene <- sig_gene_df |>
  pull(ensembl)

up_gene <- sig_gene_df |>
  filter(
    logFC_scz > 0
  ) |>
  pull(ensembl)

down_gene <- sig_gene_df |>
  filter(
    logFC_scz < 0
  ) |>
  pull(ensembl)

# Export the genes for website
# query_string <- paste0(sig_gene, collapse = "::")


# Up-reg genes -----
up_ego <- enrichGO(
  gene = up_gene,
universe = gene_df$ensembl,
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
      "processed-data/rds/10_dx_deg_adjust_spd",
      "GO_up_gene.csv")
    )

# Down-reg genes ----
down_ego <- enrichGO(
  gene = down_gene,
universe = gene_df$ensembl,
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
      "processed-data/rds/10_dx_deg_adjust_spd",
      "GO_down_gene.csv")
    )



# GO enrichment (Overrepresentative analysis) ----
gc_list <- list(
  "All" = sig_gene,
  "Up" = up_gene,
  "Down" = down_gene
)

GO_ora <- compareCluster(
  geneCluster = gc_list,
  fun = enrichGO,
  universe = gene_df$ensembl,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "ENSEMBL",
  pAdjustMethod = "BH",
  # pvalueCutoff = 0.01,
  # qvalueCutoff = 0.05,
  readable = TRUE
)

dir.create(
  here("processed-data/PB_dx_genes/enrichment"),
  showWarnings = FALSE
)

GO_ora@compareClusterResult |>
  write.csv(
    here(
      "processed-data/PB_dx_genes/enrichment",
      paste0("test_GO_ora_", adj_p_cutoff, ".csv")
    )
  )


# Make dot plot ----

go_ora_list <- c("BP", "MF", "CC") |>
  set_names() |>
  imap(
    ~ compareCluster(
      geneCluster = gc_list,
      fun = enrichGO,
      universe = gene_df$ensembl,
      OrgDb = org.Hs.eg.db,
      ont = .x,
      keyType = "ENSEMBL",
      pAdjustMethod = "BH",
      # pvalueCutoff = 0.01,
      # qvalueCutoff = 0.05,
      readable = TRUE
    )
  )

pdf(
  here(
    "plots/PB_dx_genes/enrichment",
    paste0("test_GO_ora_", adj_p_cutoff, ".pdf")
  ),
  height = 10, width = 10
)
go_ora_list |> iwalk(
  ~ dotplot(
    .x,
    showCategory = 10, # Show top 10 genes of each category
    split = "NULL",
    includeAll = FALSE,
    group = FALSE,
    title = .y
  )
) + scale_x_discrete(limits = c("Down\n(92)", "Up\n((66)", "All\n((158)"))
dev.off()


# Visualization ----
tmp <- dotplot(
  GO_ora,
  showCategory = 10, # Show top 10 genes of each category
  split = "NULL",
  includeAll = FALSE,
  group = TRUE
) +
  scale_x_discrete(limits = c("Up", "Down", "All"))

## Debug why there's a category that is NA
tmp$data |> View()

# Results ----
GO_ora@compareClusterResult |> View()

GO_ora@compareClusterResult |>
  group_by(Cluster) |>
  slice_head(n = 10) |>
  View()



# Deprecated Code ----
## All genes ----
ego <- enrichGO(
  gene = down_gene,
universe = gene_df$ensembl,
OrgDb = org.Hs.eg.db,
ont = "ALL",
keyType = "ENSEMBL",
pAdjustMethod = "BH",
pvalueCutoff = 0.01,
qvalueCutoff = 0.05,
readable = TRUE
)

ego@result |> View()

edox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, node_label="all", 
        cex_label_category = 1.2) 

library(enrichplot)
edo <- pairwise_termsim(ego)
emapplot(edo)
treeplot(edo)
# # dotplot(ego)
# write.csv(
#   ego@result,
#   here(
#     "processed-data/PB_dx_genes/enrichment",
#     "GO_ERA_PRECAST_07.csv"
#   )
# )
## Up reg genes ----
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
## Down reg genes ----
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
