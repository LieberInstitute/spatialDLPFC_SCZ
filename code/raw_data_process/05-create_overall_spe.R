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

# Confirm if the number of samples match with meta
stopifnot( nrow(expr_meta) == length(unique(spe$sample_id)))

spe$ManualAnnotation <- NULL

# Only for testing purpose
# spe <- readRDS(here("processed-data/rds/spe","01_build_spe/", "spe_raw.rds"))

## TODO: Add the experimental information
### TODO: read in meta information
source(here("code", "raw_data_process", "import_dx.R"))
# Save the dx data as the meta data
metadata(spe) <- clean_df


# col_df <- colData(spe) |> data.frame() |> 
#   mutate(brain_num = str_remove(sample_id , "_[a-zA-Z][0-9]")) |>
#   left_join(clean_df, by = "brain_num")


# spe$key <- paste0(colnames(spe), "_", spe$sample_id) # In spatialLIBD::read10xVisiumWrapper
# spe$subject <- sample_info$subjects[match(spe$sample_id, sample_info$sample_id)]
# spe$region <- sample_info$regions[match(spe$sample_id, sample_info$sample_id)]
# spe$sex <- sample_info$sex[match(spe$sample_id, sample_info$sample_id)]
# spe$age <- sample_info$age[match(spe$sample_id, sample_info$sample_id)]
# spe$diagnosis <- sample_info$diagnosis[match(spe$sample_id, sample_info$sample_id)]
# spe$sample_id_complete <- spe$sample_id
# spe$sample_id <- gsub("_2", "", spe$sample_id)




## Add information used by spatialLIBD
# Seems incorprated in spatialLIBD::read10xVisiumWrapper
# is_mito <- which(seqnames(spe) == "chrM")
# spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
# spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi

# TODO: add this information
## Read in cell counts and segmentation results
segmentations_list <- lapply(#unique(spe$sample_id), 
  unique(expr_meta$sample_name), # TEST only
                             function(sampleid) {
  
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
## Merge them (once the these files are done, this could be replaced by an rbind)
# segmentations <- Reduce(function(...) merge(..., all = TRUE), segmentations_list[lengths(segmentations_list) > 0])

# ## Add the information
# segmentation_match <- match(spe$key, segmentations$key)
# segmentation_info <- segmentations[segmentation_match, -which(colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")), drop = FALSE]
# colData(spe) <- cbind(colData(spe), segmentation_info)

# vis_grid_gene(
#   spe = spe,
#   geneid = "Ndata",
#   pdf = here::here("plots", "01_build_spe", "seg_count.pdf"),
#   assayname = "counts",
#   auto_crop = FALSE,
#   spatial = FALSE
# )

# mean(spe$count)
# pdf(here::here("plots", "01_build_spe", "cells_per_spot.pdf"))
# boxplot(spe$count)
# dev.off()

## Remove genes with no data
# no_expr <- which(rowSums(counts(spe)) == 0)

# length(no_expr)
# [1] 6936
# length(no_expr) / nrow(spe) * 100
# [1] 18.9503


# spe_raw <- spe

dir.create(
  here::here("processed-data", "rds", 
             "spe", "01_build_spe"),
  recursive = T
)

saveRDS(spe_raw, 
        file = here::here("processed-data", "rds", 
                          "spe", "01_build_spe", "spe_raw.rds"))




