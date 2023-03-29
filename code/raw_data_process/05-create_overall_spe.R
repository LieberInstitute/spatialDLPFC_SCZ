library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)


# TODO: replace this part with the preprocess_meta_file

sample_info <- data.frame(
  sample_id = list.files(
    here("processed-data", "spaceranger"),
    pattern = "^Br"
  )
) |> 
  mutate(
    sample_path = file.path(
      here::here("processed-data", "spaceranger"),
      sample_id,
      "outs"
    )
  )

stopifnot(all(file.exists(sample_info$sample_path)))

## Related to https://github.com/drighelli/SpatialExperiment/issues/135
stopifnot(packageVersion("SpatialExperiment") >= "1.9.5")
## Related to https://github.com/LieberInstitute/spatialLIBD/commit/1ff74a8be040bf4cdbc906700913ad12fc929232
stopifnot(packageVersion("spatialLIBD") >= "1.11.10")

## Build SPE object
Sys.time()
spe <- spatialLIBD::read10xVisiumWrapper(
  sample_info$sample_path,
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = c("lowres", "hires", "detected", "aligned"),
  load = TRUE,
  reference_gtf = "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
)
Sys.time()

# TODO: Confirm if the number of samples match with meta


## TODO Add the experimental information
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
# ## Read in cell counts and segmentation results
# segmentations_list <- lapply(sample_info$sample_id, function(sampleid) {
#   current <- sample_info$sample_path[sample_info$sample_id == sampleid]
#   file <- file.path(current, "spatial", "tissue_spot_counts.csv")
#   if (!file.exists(file)) {
#     return(NULL)
#   }
#   x <- read.csv(file)
#   x$key <- paste0(x$barcode, "_", sampleid)
#   return(x)
# })
# ## Merge them (once the these files are done, this could be replaced by an rbind)
# segmentations <- Reduce(function(...) merge(..., all = TRUE), segmentations_list[lengths(segmentations_list) > 0])
# 
# ## Add the information
# segmentation_match <- match(spe$key, segmentations$key)
# segmentation_info <- segmentations[segmentation_match, -which(colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")), drop = FALSE]
# colData(spe) <- cbind(colData(spe), segmentation_info)

vis_grid_gene(
  spe = spe,
  geneid = "Ndata",
  pdf = here::here("plots", "01_build_spe", "seg_count.pdf"),
  assayname = "counts",
  auto_crop = FALSE,
  spatial = FALSE
)

# mean(spe$count)
# pdf(here::here("plots", "01_build_spe", "cells_per_spot.pdf"))
# boxplot(spe$count)
# dev.off()

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)

length(no_expr)
# [1] 6936
length(no_expr) / nrow(spe) * 100
# [1] 18.9503

spe_raw <- spe

dir.create(
  here::here("processed-data", "rds", 
                      "spe", "01_build_spe"),
           recursive = T
  )

saveRDS(spe_raw, 
     file = here::here("processed-data", "rds", 
                       "spe", "01_build_spe", "spe_raw.rds"))



## Remove spots without counts
spe <- spe_raw[-no_expr, ]
dim(spe)
# [1] 29665 39936

pryr::object_size(spe)
# 1.81 GB

spe <- spe[, which(colData(spe)$in_tissue)]
dim(spe)
# [1] 29665 30858


## Remove spots without counts
# NOTE: Compare to spatilDLPFC
# all in$tissue spots have reads
# n_spot_no_counts <- colSums(counts(spe)) == 0
spe <- spe[, colSums(counts(spe)) != 0]
dim(spe)
# [1] 29665 30858


dir.create(
  here::here("plots", "01_build_spe"),
  recursive = T
)

## Inspect in vs outside of tissue
vis_grid_clus(
  spe = spe_raw,
  clustervar = "in_tissue",
  pdf = here::here("plots", "01_build_spe", "in_tissue_grid.pdf"),
  sort_clust = FALSE,
  colors = c("TRUE" = "grey90", "FALSE" = "orange")
)

## Roughly speaking, looks very nice.


# TODO: how many proportion of the spots have missing

# summary(spe_raw$sum_umi[which(!colData(spe_raw)$in_tissue)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0    66.0   170.0   291.6   326.0 39053.0

# hist(log10(spe_raw$sum_umi[which(!colData(spe_raw)$in_tissue)]))
# 0  1  2  3  4  5
# 12 17 50 59 65 81


# TODO: error
vis_grid_gene(
  spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
  geneid = "sum_umi",
  pdf = here::here("plots", "01_build_spe", "out_tissue_sum_umi_all.pdf"),
  assayname = "counts",
  auto_crop = FALSE
)
# Error in frame_lims$y_min:frame_lims$y_max : NA/NaN argument
# Auto_crop doesn't work in this case.

summary(spe_raw$expr_chrM_ratio[which(!colData(spe_raw)$in_tissue)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#  0.0000  0.1500  0.1947  0.2072  0.2500  1.0000      12

vis_grid_gene(
  spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "01_build_spe", "out_tissue_expr_chrM_ratio_all.pdf"),
  assayname = "counts",
  auto_crop = FALSE
)

vis_grid_gene(
  spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
  geneid = "sum_gene",
  pdf = here::here("plots", "01_build_spe", "out_tissue_sum_gene_all.pdf"),
  assayname = "counts"
)

summary(spe$sum_umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1    1358    2310    2674    3540   46881
vis_grid_gene(
  spe = spe,
  geneid = "sum_umi",
  pdf = here::here("plots", "01_build_spe", "in_tissue_sum_umi_all.pdf"),
  assayname = "counts"
)

summary(spe$sum_gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1     910    1426    1508    1999    8344

vis_grid_gene(
  spe = spe,
  geneid = "sum_gene",
  pdf = here::here("plots", "01_build_spe", "in_tissue_sum_gene_all.pdf"),
  assayname = "counts"
)

summary(spe$expr_chrM_ratio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.07279 0.10440 0.12106 0.15513 1.00000
vis_grid_gene(
  spe = spe,
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "01_build_spe", "in_tissue_expr_chrM_ratio_all.pdf"),
  assayname = "counts"
)

vis_grid_gene(
  spe = spe_raw,
  geneid = "sum_umi",
  pdf = here::here("plots", "01_build_spe", "all_sum_umi.pdf"),
  assayname = "counts"
)

vis_grid_gene(
  spe = spe_raw,
  geneid = "sum_gene",
  pdf = here::here("plots", "01_build_spe", "all_sum_gene.pdf"),
  assayname = "counts"
)

vis_grid_gene(
  spe = spe_raw,
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "01_build_spe", "all_expr_chrM_ratio.pdf"),
  assayname = "counts"
)

