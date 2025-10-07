# Load libraries ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
  library(ggridges)
  library(ggpubr)
})

# Load data ----
## Load layer-restricted DEGs -----
spd_deg_df <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_restricted_degs_all_spds.csv"
  )
)


## Load Ruzicka DEGs -----
ruzicka_deg_df <- readRDS(
  # NOTE: this is a list of cell type specific data frames
  file = here(
    "processed-data/rds/12_cross-study_enrichment",
    "ruzicka_cell_type_DEG.rds"
  )
) |> bind_rows()


# Visualization ----
## layer-strictited DEGs ----


## ridge plot
# Order the y-axis by PRECAST_spd factor levels (e.g., by frequency or custom order)
spd_deg_df <- spd_deg_df %>%
  mutate(PRECAST_spd = factor(PRECAST_spd,
    levels = c(
      "SpD07-L1",
      "SpD06-L2/3",
      "SpD02-L3/4",
      "SpD05-L5",
      "SpD03-L6",
      "SpD01-WMtz",
      "SpD04-WM"
    ) |> rev(),
    labels = c(
      "SpD07-L1/M",
      "SpD06-L2/3",
      "SpD02-L3/4",
      "SpD05-L5",
      "SpD03-L6",
      "SpD01-WMtz",
      "SpD04-WM"
    ) |> rev()
  ))

p_layer_restr <- ggplot(
  spd_deg_df,
  aes(
    x = abs(logFC), y = PRECAST_spd, fill = PRECAST_spd
  )
) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, color = "white") +
  labs(
    title = "Layer-restricted DEGs",
    x = "|log2(FC in SCZ)|",
    y = "PRECAST\n Spatial Domain"
  ) +
  scale_fill_manual(
    values = setNames(
      Polychrome::palette36.colors(13)[seq.int(7)],
      c(
        "SpD07-L1",
        "SpD06-L2/3",
        "SpD02-L3/4",
        "SpD05-L5",
        "SpD03-L6",
        "SpD01-WMtz",
        "SpD04-WM"
      )
    )
  ) +
  theme_classic() +
  theme(legend.position = "none")


## Ruzicka DEGs ----
p_ruzicka_DE <- ggplot(
  ruzicka_deg_df,
  aes(
    x = abs(Meta_logFC), y = cell_type, fill = cell_type
  )
) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, color = "white") +
  labs(
    title = "Ruzicka et al.",
    x = "|log2(FC in SCZ)|\n (truncated at 0.5)",
    y = "Cell Type"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 0.5))


## Assemble panneled plot ----

ggarrange(
  p_layer_restr,
  p_ruzicka_DE,
  ncol = 2,
  labels = c("A", "B")
) |>
  ggsave(
    filename = here(
      "plots/discussion",
      "ridge_plot_effect_size_comparison.pdf"
    )
  )

ggarrange(
  p_layer_restr,
  p_ruzicka_DE,
  ncol = 2,
  labels = c("A", "B")
) |>
  ggsave(
    filename = here(
      "plots/discussion",
      "ridge_plot_effect_size_comparison.png"
    )
  )

# Session info ----
session_info()
