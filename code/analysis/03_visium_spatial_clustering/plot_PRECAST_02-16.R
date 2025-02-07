# Load library -----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(escheR)
  library(sessioninfo)
})

# Load data ----
## SPE data ----
spe <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

# error prevention
stopifnot(all(spe$in_tissue))

# create sample_label
spe$sample_label <- paste0(
  spe$brnum, "_", toupper(spe$dx)
)

## Load PRECAST df ----
PRECAST_df <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

## Merge PRECAST df ----
precast_vars <- grep(
  "^PRECAST_", colnames(PRECAST_df),
  value = TRUE
)

spe <- spe[, spe$key %in% PRECAST_df$key]
col_data_df <- PRECAST_df |>
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )
rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

# error prevention
stopifnot(is.character(spe$PRECAST_07))

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

# Make plots ----
## Iterat over all possible PRECAST settings
# More data driven approach
# colnames(colData(spe)) |>
#   grep(pattern = "PRECAST_", x = _, value = TRUE) |>

panel_list <-
  set_names(2:13, ~ sprintf("PRECAST_%02s", .)) |>
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
              Polychrome::palette36.colors(13)[seq.int(.k)],
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
  nrow = 6, ncol = 2 # ,
  # labels = sprintf("k=%02d", 2:13)
  # TODO: use the max precast label
  # common.legend = TRUE,
  # legend = "bottom"
) |> ggsave(
  here(
    "plots/03_visium_spatial_clustering",
    "spot_plot_PRECAST_2-13_rep_samples.pdf"
  ),
  plot = _,
  height = 11, width = 7
)



# Session Info ----
sessioninfo()
