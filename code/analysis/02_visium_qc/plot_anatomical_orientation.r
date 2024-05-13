# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(spatialLIBD)
  library(here)
  library(scater)
  library(escheR)
  library(sessioninfo)
})


# Load Data ----
## Load Spe ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "raw_spe_wo_SPG_N63.rds"
  )
)

## Remove out_tissue_spots ---
spe <- spe[, spe$in_tissue == TRUE]

## Calcualte normaled counts ----
spe <- spe[, colSums(counts(spe)) != 0]
spe <- logNormCounts(spe)

## Adjust sample name to Brnum_dx ---
spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$brnum <- metadata(spe)$dx_df$subject[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$sample_id <- paste0(
  spe$brnum, "_", spe$dx
)


# Create visualization ----

plot_list <- vector("list", length = 63) |>
  setNames(
    unique(spe$sample_id)
  )

for (.smp in names(plot_list)) {
  # browser()
  sub_spe <- spe[, spe$sample_id == .smp]

  RBG_norm_MBP <- scales::rescale(
    # TODO: change this to logcounts
    logcounts(sub_spe)[which(rowData(sub_spe)$gene_name == "MBP"), ] |>
      scale(center = TRUE, scale = FALSE) |>
      sapply(max, 0),
    to = c(0, 1)
  )
  RBG_norm_PCP4 <- scales::rescale(
    logcounts(sub_spe)[which(rowData(sub_spe)$gene_name == "PCP4"), ] |>
      scale(center = TRUE, scale = FALSE) |>
      sapply(max, 0),
    to = c(0, 1)
  )
  RBG_norm_SNAP25 <- scales::rescale(
    logcounts(sub_spe)[which(rowData(sub_spe)$gene_name == "SNAP25"), ] |>
      scale(center = TRUE, scale = FALSE) |>
      sapply(max, 0),
    to = c(0, 1)
  )

  sub_spe$rbg_val <- rgb(
    red = RBG_norm_MBP,
    green = RBG_norm_PCP4,
    blue = RBG_norm_SNAP25
  )

  plot_list[[.smp]] <- make_escheR(sub_spe) |>
    add_fill("rbg_val") +
    scale_fill_identity() +
    theme(
      plot.background = element_rect(fill = "black"),
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 40)
    )
}




# spe$rbg_val <- rgb(
#   red = RBG_norm_MBP,
#   green = RBG_norm_PCP4,
#   blue = RBG_norm_SNAP25
# )

# plot_list <- vis_grid_gene(
#   spe,
#   gene_id = "rbg_val",
#   sample_order = unique(spe$sample_id) |> sort(),
#   spatial = FALSE,
#   point_size = 2,
#   return_plots = TRUE
# )

pdf(
  file = here::here(
    "plots", "02_visium_qc",
    paste0("rgb_plot_anatomy_sfigur.pdf")
  ),
  height = 6 * 8, width = 6 * 8 - 1
  # height = 1 * 8, width = 6 * 8 - 1 # Test
)
# Plot png for place holder in google drive.
# png(
#   file = here::here(
#     "plots", "02_visium_qc",
#     paste0("vis_clus_sample_aware_low_lib_size_sfigur.png")
#   ),
#   height = 5 * 8, width = 8 * 8
# )

# Plot ntc and scz in separate groups
idx_list <- spe$dx |>
  unique() |>
  set_names() |>
  imap(~ str_detect(names(plot_list), .x))

for (.idx in idx_list) {
  print(
    ggpubr::ggarrange(
      plotlist = plot_list[.idx],
      ncol = 6, nrow = 6,
      common.legend = TRUE
    )
  )
}

# Test code (TODO: delete)
# ggpubr::ggarrange(
#   plotlist = plot_list,
#   ncol = 6, nrow = 1,
#   common.legend = TRUE
# )

dev.off()

# Session Info
session_info()
