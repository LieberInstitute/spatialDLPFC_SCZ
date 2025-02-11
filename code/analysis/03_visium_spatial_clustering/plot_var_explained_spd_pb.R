# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
})

# spd_rds <- list.files(
#   here(
#     "processed-data", "rds", "layer_spd"
#   ),
#   pattern = ".rds"
# )

# PCA analysis ----

# dir.create(
#   path = here(
#     "plots/PB_DE", "SpD_PB_PCA"
#   ),
#   showWarnings = FALSE
# )

# for (.file in spd_rds) {
# .file <- spd_rds[1]

# .spd <- str_remove(.file, "test_spe_pseudo_") |>
#   str_remove(".rds")
.spd <- "PRECAST_07"

sce_pseudo <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    paste0("test_spe_pseudo_", .spd, ".rds")
  )
)


set.seed(20240411)
sce_pseudo <- runPCA(sce_pseudo)

pdf(
  here(
    "plots/PB_DE", "SpD_PB_PCA",
    paste0("test_", .spd, ".pdf")
  )
)
# sce_pseudo <- runPCA(sce_pseudo)
# for (.var in c("age", "sex", "spd", "dx", "slide_id", "lot_num")) {
#   plotPCA(
#     sce_pseudo,
#     colour_by = .var,
#     ncomponents = 6,
#     point_size = 0.3,
#     label_format = c("%s %02i", " (%i%%)")
#   ) |> print()
# }


dat_mat <- getVarianceExplained(sce_pseudo,
  variables = c(
    "dx", "age", "sex", "slide_id", "PRECAST_07" # , "sample_id")
  )
)

pdf(
  here("plots/03_visium_spatial_clustering",
  "den_plot_var_explained_by_spd.pdf"),
  height = 5
)
plotExplanatoryVariables(sce_pseudo,
  variables = c(
    "dx", "age", "sex", "slide_id", "PRECAST_07" # , "sample_id")
  )
)
dev.off()
# }


# Session Info ----
sessioninfo::session_info()
