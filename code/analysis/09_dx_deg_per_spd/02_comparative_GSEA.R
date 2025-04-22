# Load required packages ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  # library(DOSE)
})


# Load SpD stratified DEGs ----
# create a list of spd-stratified csv files
spd_files <- list.files(
  here("processed-data/rds/09_dx_deg_per_spd"),
  pattern = "dx-deg_PRECAST07_.*\\.csv",
  full.names = TRUE
)

spd_deg_list <- map(
  spd_files,
  ~ read_csv(.x, col_types = cols())
)

names(spd_deg_list) <- str_extract(
  spd_files,
  "(?<=dx-deg_PRECAST07_).*?(?=\\.csv)"
)

## Understand the number of sig gnenes ----
# At the nominal p-value cutoff of 0.05
map_int(
  spd_deg_list,
  ~ sum(.x$p_value_scz <= 0.05)
)
# SpD01_WMtz SpD02_L3_4   SpD03_L6   SpD04_WM   SpD05_L5 SpD06_L2_3   SpD07_L1
#        396        975        636        671        932        794        668

# At the adjusted p-value cutoff of 0.10
map_int(
  spd_deg_list,
  ~ sum(.x$fdr_scz <= 0.10)
)
# SpD01_WMtz SpD02_L3_4   SpD03_L6   SpD04_WM   SpD05_L5 SpD06_L2_3   SpD07_L1
#          0          0          0          0          0          0          1

# Run clusterProfiler GSEA ----
# NOTE: https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html

map_int(
  spd_deg_list, ~ nrow(.x)
)

# Toy example ----
spd_gene_list <- map(spd_deg_list, ~ {
  .x |>
    filter(.x$p_value_scz <= 0.05) |>
    pull(ensembl)
})

GO_ora <- compareCluster(
  geneCluster = spd_gene_list,
  fun = "enrichGO",
  universe = spd_deg_list[[1]]$ensembl,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENSEMBL",
  pAdjustMethod = "BH",
  # pvalueCutoff = 0.01,
  # qvalueCutoff = 0.05,
  readable = TRUE
)


write_csv(
  GO_ora@compareClusterResult,
  here("code/analysis/09_dx_deg_per_spd/tmp_GO_ora_BP_only.csv")
)


# gsea analysis ----
tmp_res <- spd_deg_list |>
  imap(.f = function(.df, idx) {

    # format gene list
    geneList <- .df |>
      arrange(desc(t_stat_scz)) |>
      select(ensembl, t_stat_scz) |>
      deframe()

    set.seed(20250422)
    # Run GSEA
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

    # make gene readable
    ego_gsea@result <- ego_gsea@result |>
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

    # save results
    write_csv(
      ego_gsea@result,
      here(
        "code/analysis/09_dx_deg_per_spd",
        paste0("gsea_PRECAST07_", gsub("[^[:alnum:]_]", "_", idx), ".csv")
      )
    )
  })



dotplot(GO_ora)

cnetplot(GO_ora)


str(GO_ora)
GO_ora@compareClusterResult |> View()


# Make the comparative plot


# Session Info ----
sessioninfo::session_info()
