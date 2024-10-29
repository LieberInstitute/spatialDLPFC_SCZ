# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(spatialLIBD)
  library(limma)
  library(sessioninfo)
  library(here)
  library(glue)
})

# Loop over all SPG + analysis level ----
analysis_type <- c("SPD", "donor")
# spg_names <- c("pnn_pos", "neuropil_pos", "neun_pos", "vasc_pos")
spg_names <- c("pnn_N_neighbors")
# spg_names <- paste0("pnn_pos_", c(#"L2345", #"L23", 
# "L34", "L5"))

analysis_combo <- expand.grid(
  type = analysis_type,
  spg = spg_names
)

# NOTE:
# only use when testing the manually created SPD subsetted data
# `pnn_L2345\viz_spd_pb_data.R`
# analysis_combo <- tibble(
#   type = c("SPD", "donor", "donor", "donor"),
#   spg = c("pnn_pos_L2345", "pnn_pos_L23", "pnn_pos_L34", "pnn_pos_L5")
# )



pb_files <- analysis_combo |> glue_data("test_{type}_pseudo_{spg}.rds")

pb_files |>
  walk(.f = function(.pb_file) {
    # Error prevention
    # Check if pseudobulk rds file exists
    stopifnot(
      file.exists(
        here(
          "processed-data", "rds", "PB_dx_spg",
          .pb_file
        )
      )
    )

    ## Load Data ----
    sce_pseudo <- readRDS(
      here(
        "processed-data", "rds", "PB_dx_spg",
        .pb_file
      )
    )
    # colData(sce_pseudo) |> str()
    ## Run pipeline ----

    if (grepl("donor", .pb_file)) {
      ## At the donor level
      cat("Donor level DE analysis\n")
      dx_res <- registration_stats_enrichment(
        sce_pseudo,
        block_cor = NaN,
        covars = c("age", "sex"),
        var_registration = "dx",
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
      )
    } else if (grepl("SPD", .pb_file)) {
      cat("SPD level DE analysis\n")
      ## Adjust for SPD
      dx_mod <-
        registration_model(
          sce_pseudo,
          covars = c("age", "sex"),
          var_registration = "dx"
        )

      dx_block_cor <- registration_block_cor(
        sce_pseudo,
        registration_model = dx_mod,
        var_sample_id = "sample_id"
      )

      dx_res <- registration_stats_enrichment(
        sce_pseudo,
        block_cor = dx_block_cor,
        covars = c("age", "sex", "PRECAST_07"),
        var_registration = "dx",
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
      )
    } else {
      stop("No analysis for this type of PB data.")
    }


    ### Save the DE result ----
    write_csv(
      dx_res,
      here(
        "processed-data/spg_pb_de",
        .pb_file |> str_replace(".rds", ".csv")
      )
    )
  })

# Session Info ----
session_info()
