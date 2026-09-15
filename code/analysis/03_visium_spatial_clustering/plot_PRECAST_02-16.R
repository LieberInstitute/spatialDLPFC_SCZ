# Load library -----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(escheR)
  library(ggpubr)
  library(sessioninfo)
})

# Load data ----
## SPE data ----
spe <- readRDS(
  here(
    "processed-data/rds/03_visium_spatial_clustering",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

# error prevention
stopifnot(all(spe$in_tissue))

# create sample_label
spe$sample_label <- paste0(
  spe$brnum, "_", toupper(spe$dx)
)

## Subset to representative samples ----
rep_sample_id <- c("V13M06-342_D1", "V13M06-343_D1")

# error prevention
stopifnot(
  all(
    rep_sample_id %in%
      unique(spe$sample_id)
  )
)

spe <- spe[, spe$sample_id %in% rep_sample_id]

# error prevention
stopifnot(
  length(
    unique(spe$sample_id)
  ) == 2
)

## Load and add PRECAST labels ----
PRECAST_df <- readRDS(
  here(
    "processed-data/rds/03_visium_spatial_clustering",
    "PRECAST_label_df_semi_sup_k_2-16.rds"
  )
)

precast_vars <- grep(
  "^PRECAST_", colnames(PRECAST_df),
  value = TRUE
)
precast_idx <- match(spe$key, PRECAST_df$key)

stopifnot(!anyNA(precast_idx))

colData(spe)[, precast_vars] <- DataFrame(
  PRECAST_df[precast_idx, precast_vars, drop = FALSE]
)

# error prevention
stopifnot(is.character(spe$PRECAST_07))

rm(PRECAST_df)
gc()

# Make plots ----
## Iterat over all possible PRECAST settings
# More data driven approach
# colnames(colData(spe)) |>
#   grep(pattern = "PRECAST_", x = _, value = TRUE) |>

panel_list <-
  set_names(2:16, ~ sprintf("PRECAST_%02s", .)) |>
  map(.f = function(.k) {
    # iterate over samples
    precast_k <- sprintf("PRECAST_%02d", .k)
    p_list <- rep_sample_id |>
      map(.f = function(.sample) {
        sub_spe <- spe[, spe$sample_id == .sample]

        ret_p <- make_escheR(sub_spe) |>
          # TODO: adjust point size here.
          add_fill(
            var = precast_k,
            point_size = 0.8
          ) +
          # labs(title = unique(sub_spe$sample_label)) +
          # Adjust the color palette
          scale_fill_manual(
            name = "Spatial Domain",
            values = set_names(
              Polychrome::palette36.colors(16)[seq.int(.k)],
              unique(sub_spe[[precast_k]]) |> sort()
            )
          ) +
          theme(
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            plot.margin = margin(t = 10, b = 10)
          )
        return(ret_p)
      })

    # browser()
    p_panel <- ggarrange(
      plotlist = p_list,
      nrow = 1,
      ncol = 2,
      common.legend = TRUE, legend = "none"
    ) #|>
    # annotate_figure(
    #   fig.lab = sprintf("k=%02d", .k),
    #   fig.lab.size = 12,
    #   fig.lab.pos = "top"
    # )
    return(p_panel)
  })

## Format individual plots ----

# Make panels ----
cowplot::plot_grid(
  plotlist = panel_list,
  nrow = 8, ncol = 2 # ,
  # labels = sprintf("k=%02d", 2:16)
  # TODO: use the max precast label
  # common.legend = TRUE,
  # legend = "bottom"
) |> ggsave(
  here(
    "plots/03_visium_spatial_clustering",
    "spot_plot_PRECAST_2-16_rep_samples.pdf"
  ),
  plot = _,
  height = 14.7, width = 7
)



# Session Info ----
session_info()
