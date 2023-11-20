# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(tidyverse)
library(sessioninfo)

## Related to 
## https://github.com/drighelli/SpatialExperiment/issues/135
stopifnot(packageVersion("SpatialExperiment") >= "1.9.5")
## Related to 
## https://github.com/LieberInstitute/spatialLIBD/commit/1ff74a8be040bf4cdbc906700913ad12fc929232
stopifnot(packageVersion("spatialLIBD") >= "1.11.10")



# Load SPE Object ---------------------------------------------------------
# TODO: add path
path_raw_spe <- "REPLACE_PATH HERE"

# TODO: check if spe exists
stopifnot(
  "Please create raw spe object first" =
    file.exists(path_raw_spe)
)

# TODO: load in spe object
raw_spe <- readRDS(path_raw_spe)

# TODO: error prevention
stopifnot(
  "Non DLPFC sample still exists" = 
    length(
      raw_spe$sample_id |> unique()
    ) == 63
)

# Load SPG file path ---------------------------------------------
raw_expr_meta <- read.csv(
  here("code", "raw_data_process",
       "sample_meta_path.csv"),
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
  filter(sample_name !=  "V12F14-053_B1")

stopifnot(
  "Non DLPFC sample still exists" = 
    length(
      expr_meta$sample_name |> unique()
    ) == 63
)

# Add SPG Channels --------------------------------------------------------
## Read in cell counts and segmentation results (TODO:)
spg_df <- map_dfr(
  unique(raw_spe$sample_id),
  function(sampleid) {
    # browser()
    current <- expr_meta$sr_fldr_path[expr_meta$sample_name == sampleid]
    stopifnot("more than one match or no match" = length(current)==1)
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
  })

stopifnot("SPG masking data doesn't match spe spots 1-to-1" = nrow(spg_df) == ncol(raw_spe))

colnames(spg_df)
# [1] "barcode"    "tissue"     "row"        "col"        "imagerow"
# [6] "imagecol"   "NAF"        "PAF"        "CNAF"       "NClaudin5"
# [11] "PClaudin5"  "CNClaudin5" "NDAPI"      "PDAPI"      "CNDAPI"
# [16] "NNeuN"      "PNeuN"      "CNNeuN"     "NWFA"       "PWFA"
# [21] "CNWFA"      "key"

# TODO: explain
# Note: 4 channels:
#* AF -
#* Claudin5 -
#* Neun -
#* WFA -
# Note 2: Naming Convention
#* N-channel-:
#* P-channel-:
#* CN-channel-:


# TODO: clean up names for seg_df
# name starts with SPG_
# ends with P and N
colnames(spg_df) <- paste0("spg_", colnames(spg_df))

spe <- raw_spe

col_data_df <- colData(spe) |> data.frame() |>
  left_join(
    spg_df, by = c("key" = "spg_key"),
    relationship = "one-to-one"
  ) |>
  left_join(
    expr_meta |>
      select(
        sample_name,
        BrNumbr),
    by = c("sample_id" = "sample_name")
  )

# Add the information
colData(spe) <- DataFrame(col_data_df)  # Will remove colnames(spe)
colnames(spe) <- spe$key
stopifnot(any(is.na(spe$spg_CNClaudin5)))

# * Save SPE-SPG object ---------------------------------------------------

stopifnot(!is.null(colnames(spe)))

saveRDS(
  spe,
  file = here::here("processed-data", "rds",
                    "spe", "01_build_spe",
                    # TODO: rename
                    "test_raw_spe_spg_36.rds")
)


# Session Info ------------------------------------------------------------
session_info()