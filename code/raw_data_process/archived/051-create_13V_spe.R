library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)



# Read In Data ------------------------------------------------------------

expr_meta <- readr::read_csv(
  here("code", "raw_data_process",
       "sample_meta_path.csv")
) |> filter(`Sample #` == "13v")

stopifnot(
  nrow(expr_meta) == 1
)

## Related to https://github.com/drighelli/SpatialExperiment/issues/135
stopifnot(packageVersion("SpatialExperiment") >= "1.9.5")
## Related to https://github.com/LieberInstitute/spatialLIBD/commit/1ff74a8be040bf4cdbc906700913ad12fc929232
stopifnot(packageVersion("spatialLIBD") >= "1.11.10")

## Build SPE object
Sys.time()
expr_meta <- expr_meta |> filter(`Sample #` == "13v")
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
stopifnot(ncol(spe) == length(unique(spe$sample_id))*4992)
stopifnot(nrow(expr_meta) == length(unique(spe$sample_id)))
spe$ManualAnnotation <- NULL

# Set names for each spot
spe$key <- paste0(colnames(spe), "_", spe$sample_id) # In spatialLIBD::read10xVisiumWrapper

# Add log normalization -----
library(scater)
# Create logcounts
spe <- logNormCounts(spe)

# Save Object -----
# dir.create(
#   here::here("processed-data", "rds", 
#              "spe", "01_build_spe"),
#   recursive = T, showWarnings = FALSE
# )
# 
# 
# stopifnot(!is.null(colnames(spe)))
# 
# saveRDS(spe,
# file = here::here("processed-data", "rds",
#                   "spe", "01_build_spe", "spe_raw_13v.rds"))


fldr_tissue_plots <- here("plots", "02_visium_qc",
                          "in_tissue_plots")

# Visualizations ----------------------------------------------------------


# *In-tissue spot annotation ---------------------------------------------
sample_name <- expr_meta$sample_name
vis_clus(
  spe,
  clustervar = "in_tissue",
  colors = c(
    "TRUE" = "transparent",
    "FALSE" = "grey90"#, "TRUE" = "orange"
  ) 
) |> 
  ggsave(
    filename = file.path(
      fldr_tissue_plots,
      paste0(sample_name,"_intissue",".pdf")
    ),
    height = 8,
    width = 8
  )

# * Counts ---------------------------------------------
library(escheR)

vis_gene(
  spe = spe[,spe$in_tissue == TRUE],
  geneid = "sum_umi"
)

spe$sum_log_umi <- colSums(logcounts(spe))

make_escheR(spe[, spe$in_tissue == TRUE]) |> 
  add_fill(var = "sum_log_umi")


make_escheR(spe[, spe$in_tissue == TRUE]) |> 
  add_fill(var = "sum_gene")


