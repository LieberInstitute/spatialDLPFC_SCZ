# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
})


# Load Data ----
## Load PRECAST SPE----
raw_spe <- readRDS(
  here(
    "processed-data/rds/01_build_spe",
    "raw_spe_wo_SPG_N63.rds"
  )
)

## Load PRECAST df ----
PRECAST_df <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

## Merge PRECAST df ----
precast_vars <- grep(
  "^PRECAST_", colnames(PRECAST_df),
  value = TRUE
)
raw_spe <- raw_spe[, raw_spe$key %in% PRECAST_df$key]
# raw_spe[, precast_vars] <- PRECAST_df[raw_spe$key, precast_vars]
col_data_df <- PRECAST_df |>
  right_join(
    colData(raw_spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )
rownames(col_data_df) <- colnames(raw_spe)
colData(raw_spe) <- DataFrame(col_data_df)
raw_spe <- raw_spe[, raw_spe$sample_id %in% c("V13M06-342_D1", "V13M06-343_D1")]
stopifnot(ncol(raw_spe) != 0)
gc()

## Layer Enrichment results -----
layer_rds <- list.files(
  here(
    "processed-data", "rds", "layer_enrich_test"
  ),
  pattern = ".rds"
)

## Registration to Manual Annotationn ----
manual_modeling_results <- fetch_data(type = "modeling_results")



# Layer Enrichment analysis ----
max_k <- 16

# .file <- layer_rds[[3]]
for (.file in layer_rds[1]) {
  .spd <- paste0(
    "PRECAST_",
    str_extract(.file, "(?<=PRECAST_)\\d{2}")
  )

  layer_res <- readRDS(
    here(
      "processed-data", "rds", "layer_enrich_test",
      .file
    )
  )

  ## Format enrichment test res ----
  t_stats <- layer_res[, grep("^t_stat_", colnames(layer_res))]
  colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))



  ### top 100 genes only ---
  manual_cor_g100 <- layer_stat_cor(
    t_stats,
    manual_modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = 100
  )



  annotate_registered_clusters(
    manual_cor_g100,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.1
  )



  ## Registration Heatmap ----

  ## Spatial Plots ----

  k <- unique(raw_spe[[.spd]]) |> length()

  p_list <- vis_grid_clus(
    raw_spe,
    clustervar = .spd,
    spatial = FALSE,
    sort_clust = FALSE,
    colors = set_names(
      Polychrome::palette36.colors(max_k)[seq.int(k)],
      unique(raw_spe[[.spd]]) |> sort()
    ),
    point_size = 1,
    return_plots = TRUE,
  )

  # df <- data.frame(
  #   all_spd = sprintf("SpD%02d", 1:max_k)
  # )

  # # Create a vector of 16 different colors
  # colors <- set_names(
  #   Polychrome::palette36.colors(max_k),
  #   sprintf("SpD%02d", 1:max_k)
  # )

  # # Create the plot
  # tmp_p <- ggplot(df, aes(x = all_spd, fill = all_spd)) +
  #   geom_bar() +
  #   scale_fill_manual(values = colors)

  # library(ggpubr)
  # leg <- get_legend(tmp_p)

  pdf(
    here(
      "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots",
      "reg_n_spot",
      paste0("test_reg_n_spot_", .spd, ".pdf")
    ),
    width = 10,
    height = 7
  )

  # cowplot::plot_grid(
  layer_stat_cor_plot(
    manual_cor_g100,
    max = max(manual_cor_g100)
  ) # ,
  # cowplot::plot_grid(
  cowplot::plot_grid(plotlist = p_list, ncol = 2, align = "t") |> print() # ,
  # leg,
  # ncol = 1, align = "p"
  # )
  dev.off()
}



# Session Info ----
sessioninfo::session_info()
