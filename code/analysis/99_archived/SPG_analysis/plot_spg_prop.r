# Load packages----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(escheR)
  library(sessioninfo)
})

# Load Data ----
## Load SPG spe object ----
spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_w_spg_N63.rds"
  )
)

## Load SpD data ----
finalized_spd <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

## Attach SpD label to spe ----
col_data_df <- colData(spe) |>
  data.frame() |>
  left_join(
    finalized_spd,
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

## SpD Annotation ----
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))


# Call SPG spots ----
spe$pnn_pos <- ifelse(spe$spg_PWFA > 0.05, TRUE, FALSE)
# NOTE: neuropil spot are spots doesn't have DAPI staining
spe$neuropil_pos <- ifelse(
  spe$spg_PDAPI > 0.05 & spe$spg_PDAPI < 0.5,
  FALSE, TRUE
)
spe$neun_pos <- ifelse(
  spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3,
  TRUE, FALSE
)
spe$vasc_pos <- ifelse(
  spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,
  TRUE, FALSE
)

# Viz Channels ----
spg_names <- c("pnn_pos", "neuropil_pos", "neun_pos", "vasc_pos")

## Viz Spot Plot ----

# spg_names |>
# for(.f = function(.channel) {
for (.channel in spg_names) {
  # browser()
  pdf(
    here(
      "plots/spg_analysis",
      paste0("test_spot_plot_", .channel, ".pdf")
    ),
    height = 4, width = 5
  )
  for (.sample_id in unique(spe$sample_id)) {
    p <- make_escheR(spe[, spe$sample_id == .sample_id]) |>
      add_ground("PRECAST_07") |>
      add_fill(.channel) +
      scale_color_manual(
        values = set_names(
          Polychrome::palette36.colors(7)[seq.int(7)],
          unique(spe[["PRECAST_07"]]) |> sort()
        ),
        limits = spd_anno_df$spd[order(spd_anno_df$anno_lab)],
        labels = setNames(spd_anno_df$anno_lab, spd_anno_df$spd)
      ) +
      scale_fill_manual(
        values = c(
          "TRUE" = "black",
          "FALSE" = "transparent"
        )
      ) +
      labs(title = .sample_id)
    plot(p)
  }
  dev.off()
  # }) # End of Walk
}


## Proportion of the PNN+ spots ----
col_dat <- colData(spe) |> data.frame()

pdf(
  here(
    "plots/spg_analysis",
    "test_spg_prop_per_individual.pdf"
  )
)
spg_names |>
  walk(.f = function(.channel) {
    # browser()
    p <- col_dat |>
      group_by(sample_id) |>
      summarize(
        pos_prop = sum(!!sym(.channel)) / n()
      ) |>
      left_join(
        metadata(spe)$dx_df |> select(sample_id, dx),
        by = "sample_id"
      ) |>
      ggplot() +
      geom_boxplot(
        aes(x = dx, y = pos_prop, color = dx)
      ) +
      geom_jitter(
        aes(x = dx, y = pos_prop),
        size = 0.6,
        alpha = 0.4
      ) +
      labs(title = .channel) +
      theme_minimal()

    plot(p)
  })
dev.off()


## Proportion of the PNN+ spots per layer ----

pdf(
  here(
    "plots/spg_analysis",
    "test_spg_prop_per_spd_individual.pdf"
  )
)
spg_names |>
  walk(.f = function(.channel) {
    p <- col_dat |>
      group_by(sample_id, PRECAST_07) |>
      summarize(
        prop = sum(!!sym(.channel)) / n(),
        n = n()
      ) |>
      left_join(
        metadata(spe)$dx_df |> select(sample_id, dx),
        by = "sample_id"
      ) |>
      ggplot() +
      geom_boxplot(aes(x = PRECAST_07, y = prop, color = dx)) +
      scale_x_discrete(
        limits = spd_anno_df$spd[order(spd_anno_df$anno_lab)],
        labels = setNames(spd_anno_df$anno_lab, spd_anno_df$spd)
      ) +
      labs(title = .channel) +
      theme_minimal()

    plot(p)
  })
dev.off()


# Session Info ----
session_info()
