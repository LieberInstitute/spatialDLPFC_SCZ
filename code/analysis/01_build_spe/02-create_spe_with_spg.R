# Job Script Related Notes ------------------------------------------------
# Running time: 1 hour
# Mem usage: asked for 40G, used 32GB in interactive session.

# Load Packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
})


## Related to
## https://github.com/drighelli/SpatialExperiment/issues/135
stopifnot(packageVersion("SpatialExperiment") >= "1.9.5")
## Related to
## https://github.com/LieberInstitute/spatialLIBD/commit/1ff74a8be040bf4cdbc906700913ad12fc929232
stopifnot(packageVersion("spatialLIBD") >= "1.11.10")



# Load SPE Object ---------------------------------------------------------
# Add path
path_raw_spe <- here(
  "processed-data/rds/01_build_spe",
  "raw_spe_wo_SPG_N63.rds"
)

# Check if spe exists
stopifnot(
  "Please create raw spe object first" =
    file.exists(path_raw_spe)
)

# Load in spe object
raw_spe <- readRDS(path_raw_spe)

# Error prevention
stopifnot(
  "Non DLPFC sample still exists" =
    length(
      raw_spe$sample_id |> unique()
    ) == 63
)

print("Finish loading spe")

# Load SPG file path ---------------------------------------------
raw_expr_meta <- read.csv(
  here(
    "code", "visium_data_process",
    "sample_meta_path.csv"
  ),
  header = TRUE
) |>
  filter(
    # Successful spaceranger contains outs folder
    file.exists(
      file.path(
        sr_fldr_path,
        "outs"
      )
    )
  )

expr_meta <- raw_expr_meta |>
  filter(sample_name != "V12F14-053_B1")

stopifnot(
  "Non DLPFC sample still exists" =
    length(
      expr_meta$sample_name |> unique()
    ) == 63
)

# Add SPG Channels --------------------------------------------------------
## Read in cell counts and segmentation results
spg_df <- map_dfr(
  unique(raw_spe$sample_id),
  function(sampleid) {
    current <- expr_meta$sr_fldr_path[expr_meta$sample_name == sampleid]
    stopifnot("more than one match or no match" = length(current) == 1)
    stopifnot("not the same sample" = str_detect(current, sampleid))
    # Pixel based;
    file <- file.path(current, "outs/spatial", "tissue_spot_counts.csv")
    if (!file.exists(file)) {
      warning(sampleid, " doesn't have outs/spatial/tissue_spot_counts.csv")
      return(NULL)
    }
    x <- read.csv(file)
    x$key <- paste0(x$barcode, "_", sampleid)
    return(x)
  }
)

colnames(spg_df)
# [1] "barcode"    "tissue"     "row"        "col"        "imagerow"
# [6] "imagecol"   "NDAPI"      "PDAPI"      "IDAPI"      "CNDAPI"
# [11] "NNeuN"      "PNeuN"      "INeuN"      "CNNeuN"     "NWFA"
# [16] "PWFA"       "IWFA"       "CNWFA"      "NClaudin5"  "PClaudin5"
# [21] "IClaudin5"  "CNClaudin5" "key"

# Note: 4 channels:
#* DAPI - Cells
#* Claudin5 - Blood vessels
#* Neun - Neurons
#* WFA - PNN
# Note 2: Naming Convention
#* N-channel-: Number of segmented object of interest
#* P-channel-: Percentage of pixel covered by the object of interest
#* I-channel-: Mean intensity of all the segmented pixels in the spot
#* CN-channel-: Centroid based signal count for that spot


# clean up names for seg_df
# name starts with SPG_
# ends with P and N
colnames(spg_df) <- paste0("spg_", colnames(spg_df))

cat("Finish loading spg objects\n")

spe <- raw_spe

col_data_df <- colData(spe) |>
  data.frame() |>
  left_join(
    spg_df |> select(
      -c(
        spg_tissue, spg_barcode,
        spg_row, spg_col,
        spg_imagerow, spg_imagecol
      )
    ),
    by = c(
      "key" = "spg_key"
    ),
    relationship = "one-to-one"
  ) |>
  left_join(
    expr_meta |>
      select(
        sample_name,
        BrNumbr
      ),
    by = c("sample_id" = "sample_name")
  )

# Add the information
colData(spe) <- DataFrame(col_data_df) # Remove colnames(spe)
colnames(spe) <- spe$key




# Save SPE-SPG object ---------------------------------------------------

stopifnot(!is.null(colnames(spe)))

saveRDS(
  spe,
  file = here::here(
    "processed-data", "rds", "01_build_spe",
    "raw_spe_w_spg_N63.rds"
  )
)


# Session Info ------------------------------------------------------------
session_info()
