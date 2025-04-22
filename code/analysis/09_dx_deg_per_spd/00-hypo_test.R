# Load Libray --------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(sessioninfo)
  library(tidyverse)
  library(ggrepel)
})

# Load donor-spd PB data ----
sce_pseudo <- readRDS(
  here(
    "processed-data/rds/07_dx_pseudobulk",
    "sce_pseudo_PRECAST07_donor_spd.rds"
  )
)


# Pseudo-bulk DE ----

spd_lvl <- unique(sce_pseudo$spd_label)
stopifnot(length(spd_lvl) == 7)

## Create a list of spd-stratified PB data ----
spe_list <- spd_lvl |>
  set_names() |>
  map(
    .f = function(.spd_lvl, spe) {
      # browser()
      spe[, spe$spd_label == .spd_lvl]
    },
    spe = sce_pseudo
  )

## Conduct DE per SpD ----
spe_list |>
  iwalk(
    .f = function(spe_sub, idx) {
      # spe_sub <- spd_list[[3]]
      # THis is taking the sum. The logcounts is CPM scale.
      # sce_pseudo <-
      #   registration_pseudobulk(
      #     spe_sub,
      #     var_registration = "dx",
      #     var_sample_id = "sample_id",
      #     min_ncells = 10
      #   )

      covars <- c("age", "sex")
      registration_mod <-
        registration_model(
          spe_sub,
          covars = covars,
          var_registration = "dx"
        )

      # Note: block_cor should&will return 0, as there's no block structure.
      #       It would apprear as a warning
      block_cor <-
        registration_block_cor(
          spe_sub,
          registration_model = registration_mod
        )

      results_enrichment <-
        registration_stats_enrichment(
          spe_sub,
          block_cor = block_cor,
          var_registration = "dx",
          covars = covars,
          gene_ensembl = "gene_id",
          gene_name = "gene_name"
        )

      # Save outcome ----
      write_csv(
        results_enrichment,
        here(
          "processed-data/rds/09_dx_deg_per_spd",
          # file name
          paste0(
            "dx-deg_PRECAST07_",
            # remove special characters for file name
            gsub("[^[:alnum:]_]", "_", idx),
            ".csv"
          )
        )
      )
    }
  )

# Session Info ----
sessioninfo::session_info()
