# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)
library(glue)
library(sessioninfo)


## Related to
## https://github.com/drighelli/SpatialExperiment/issues/135
stopifnot(packageVersion("SpatialExperiment") >= "1.9.5")
## Related to
## https://github.com/LieberInstitute/spatialLIBD/commit/1ff74a8be040bf4cdbc906700913ad12fc929232
stopifnot(packageVersion("spatialLIBD") >= "1.11.10")


# Confirm Space Ranger Output ---------------------------------------------
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

stopifnot(
  "Spaceranger is not successful for all 64 samples" =
    length(
      raw_expr_meta$sample_name |> unique()
    ) == 64
)


# Remove V12F14-053_B1 (Not DLPFC) ----------------------------------------
expr_meta <- raw_expr_meta |>
  filter(sample_name != "V12F14-053_B1")

stopifnot(
  "Non DLPFC sample still exists" =
    length(
      expr_meta$sample_name |> unique()
    ) == 63
)


# Build SPE Object --------------------------------------------------------
## NOTE:
## Memory space --mem=60G for 36 samples
## Required --mem=120G for 64 sample
spe <- spatialLIBD::read10xVisiumWrapper(
  samples = file.path(expr_meta$sr_fldr_path, "outs"),
  sample_id = expr_meta$sample_name,
  type = "sparse",
  data = "raw",
  images = c("lowres"), # Save memory
  # images = c("lowres", "hires", "detected", "aligned"),
  load = FALSE, # Save Memory
  # load = TRUE,
  reference_gtf = "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
)

# Confirm if the number of samples match with meta
stopifnot(nrow(expr_meta) == length(unique(spe$sample_id)))

# Confirm the number of spots are correct
if ((nrow(expr_meta) * 4992 - ncol(spe)) != 0) {
  warning(
    table(spe$sample_id) |> data.frame() |>
      filter(Freq != 4992) |>
      mutate(n_miss = 4992 - Freq) |>
      glue_data(
        "Sample {Var1} has {n_miss} missing spots."
      ) |>
      paste(sep = "\n")
  )
}
## Note:
##  Sample V13M06-281_B1 has 1 missing spots.


# Organize colData(spe) Names ---------------------------------------------
# Remove useless columns in colData(spe)
spe$ManualAnnotation <- NULL
colnames(spe) <- spe$key


# Rename Spaceranger output to sample specific (deprecated)
# colData(spe) |> names() |> grep(pattern = "10x", x = _, value = TRUE) |>
#   data.frame(name_from = _) |>
#   mutate(name_to = paste0(name_from, "_per_sample")) |>
#   # Need to do it in the global env
#   # TODO (deprecated): should be fine to use for loop
#   pmap(.f = function(name_from, name_to){
#     `$`(spe, name_to) <- `$`(spe, name_from)
#     `$`(spe, name_from) <- NULL
#   })
# (deprecated): remove reducedDim(spe)

# Add Donor Demo Info as metadata(spe) ---------
### read in meta information
source(here("code", "analysis/01_build_spe", "fun_import_dx.R"))

# Save the dx data as the meta data
metadata(spe)$dx_df <- clean_df |>
  filter(brain_num %in% expr_meta$BrNumbr) |>
  left_join(
    expr_meta |>
      select(BrNumbr, sample_name),
    by = c("brain_num" = "BrNumbr")
  ) |>
  rename(
    subject = brain_num,
    sample_id = sample_name
  ) |>
  # Remove non DLPFC sample
  dplyr::filter(
    sample_id %in% unique(spe$sample_id)
  )


# Save SPE Object --------------------------------------------------------------spe
path_spe_folder <- here::here(
  "processed-data", "rds",
  "01_build_spe"
)

dir.create(
  path_spe_folder,
  recursive = TRUE, showWarnings = FALSE
)

stopifnot(!is.null(colnames(spe)))

saveRDS(spe,
  file = here::here(
    path_spe_folder,
    "raw_spe_wo_SPG_N63.rds"
  )
)

# Session Info ------------------------------------------------------------
session_info()
