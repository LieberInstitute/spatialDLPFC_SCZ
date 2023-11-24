library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)

raw_spe <- readRDS(here("processed-data/rds/spe/01_build_spe/spe_raw.rds"))
spe <- raw_spe[, raw_spe$sample_id == "V12F14-053_A1"]

expr_meta <- read.csv(
  here("code", "raw_data_process",
       "sample_meta_path.csv"),
  header = TRUE
)

seg_df <- map_dfr(unique(spe$sample_id), 
                  function(sampleid) {
                    # browser()
                    # current <- expr_meta$sr_fldr_path[expr_meta$sample_name == sampleid]
                    # Pixel based;
                    current <- here("processed-data/spaceranger/V12F14-053_A1/")
                    # file <- file.path(current, "outs/spatial", "tissue_spot_counts.csv")
                    # centroid based,
                    file <- file.path(current, "outs/spatial", "tissue_spot_counts_centroid.csv")
                    if (!file.exists(file)) {
                      warning(sampleid, "doesn't have outs/spatial/tissue_spot_counts.csv.")
                      return(NULL)
                    }
                    x <- read.csv(file)
                    x$key <- paste0(x$barcode, "_", sampleid)
                    return(x)
                  })

stopifnot(nrow(seg_df) == ncol(spe))

colnames(seg_df) <- paste0("spg_", colnames(seg_df), "_centroid")

col_df <- colData(spe) |> data.frame() |> 
  left_join(
    seg_df, by = c("key" = "spg_key_centroid"),
    relationship = "one-to-one"
  ) 

colData(spe) <- DataFrame(col_df)

vis_gene(spe,
         geneid = "spg_NWFA_centroid",
         spatial = FALSE)
 