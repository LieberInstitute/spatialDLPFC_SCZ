library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)


# TODO: replace this part with the preprocess_meta_file

# sample_info <- data.frame(
#   sample_id = list.files(
#     here("processed-data", "spaceranger"),
#     pattern = "^V"
#   )
# ) |> 
#   mutate(
#     sample_path = file.path(
#       here::here("processed-data", "spaceranger"),
#       sample_id,
#       "outs"
#     )
#   )



expr_meta <- read.csv(
  here("code", "raw_data_process",
       "sample_meta_path.csv"),
  header = TRUE
)

stopifnot(
  all(
    file.exists(
      file.path(
        expr_meta$sr_fldr_path,
        "outs"
      )
    )
  )
)

## Related to https://github.com/drighelli/SpatialExperiment/issues/135
stopifnot(packageVersion("SpatialExperiment") >= "1.9.5")
## Related to https://github.com/LieberInstitute/spatialLIBD/commit/1ff74a8be040bf4cdbc906700913ad12fc929232
stopifnot(packageVersion("spatialLIBD") >= "1.11.10")

## Build SPE object
Sys.time()
spe <- spatialLIBD::read10xVisiumWrapper(
  samples = file.path(expr_meta$sr_fldr_path,"outs"),
  sample_id = expr_meta$sample_name,
  type = "sparse",
  data = "raw",
  images = c("lowres", "hires", "detected", "aligned"),
  load = TRUE,
  reference_gtf = "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
)
Sys.time()

# Confirm the number of spots are correct
stopifnot(ncol(spe) != length(unique(spe$sample_id))*4992)

# Confirm if the number of samples match with meta
stopifnot(nrow(expr_meta) == length(unique(spe$sample_id)))

spe$ManualAnnotation <- NULL

# Only for testing purpose
# spe <- readRDS(here("processed-data/rds/spe","01_build_spe/", "spe_raw.rds"))

# spe$key <- paste0(colnames(spe), "_", spe$sample_id) # In spatialLIBD::read10xVisiumWrapper
# spe$subject <- sample_info$subjects[match(spe$sample_id, sample_info$sample_id)]
# spe$region <- sample_info$regions[match(spe$sample_id, sample_info$sample_id)]
# spe$sex <- sample_info$sex[match(spe$sample_id, sample_info$sample_id)]
# spe$age <- sample_info$age[match(spe$sample_id, sample_info$sample_id)]
# spe$diagnosis <- sample_info$diagnosis[match(spe$sample_id, sample_info$sample_id)]
# spe$sample_id_complete <- spe$sample_id
# spe$sample_id <- gsub("_2", "", spe$sample_id)


## Read in cell counts and segmentation results
seg_df <- map_dfr(unique(spe$sample_id), 
                             function(sampleid) {
  # browser()
  current <- expr_meta$sr_fldr_path[expr_meta$sample_name == sampleid]
  file <- file.path(current, "outs/spatial", "tissue_spot_counts.csv")
  if (!file.exists(file)) {
    warning(sampleid, "doesn't have outs/spatial/tissue_spot_counts.csv.")
    return(NULL)
  }
  x <- read.csv(file)
  x$key <- paste0(x$barcode, "_", sampleid)
  return(x)
})

stopifnot(nrow(seg_df) == ncol(spe))

# TODO: clean up names for seg_df
# name starts with SPG_
# ends with P and N
colnames(seg_df) <- paste0("spg_", colnames(seg_df))

col_data_df <- colData(spe) |> data.frame() |> 
  left_join(
    seg_df, by = c("key" = "spg_key"),
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
colData(spe) <- DataFrame(col_data_df)



## Add Dx Information ---------
### read in meta information
source(here("code", "raw_data_process", "import_dx.R"))

# Save the dx data as the meta data
metadata(spe) <- clean_df |> 
  # mutate(BrNumbr = str_remove(brain_num, "Br")) |>
  filter(brain_num %in% expr_meta$BrNumbr) |> 
  # unique() |>             # TODO: to delete unique after test
  # TODO: add sample_id to the dataset
  left_join(
    expr_meta |> 
      select(BrNumbr, sample_name
            ),
    by = c("brain_num" = "BrNumbr")
  )


dir.create(
  here::here("processed-data", "rds", 
             "spe", "01_build_spe"),
  recursive = T
)

saveRDS(spe, 
        file = here::here("processed-data", "rds", 
                          "spe", "01_build_spe", "spe_raw.rds"))




