# Load library ----
suppressPackageStartupMessages({
  library(PRECAST)
  library(SpatialExperiment)
  library(scater)
  library(sessioninfo)
  library(here)
})

# Load data -----
## Load spe obejct -----
spe <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

# create sample_label
spe$sample_label <- paste0(
  spe$brnum, "_", toupper(spe$dx)
)

spe$DX <- toupper(spe$dx)

# Load PRECAST labelsdf ----
PRECAST_df <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

# Merge PRECAST df
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



## Load PRECST int object ----
seuInt <- readRDS(
  here(
    "processed-data/rds/spatial_cluster/PRECAST",
    "test_seuInt_UMAP_tsne.rds"
  )
)


## put latent embeddings  to spe object ----
# error prevention
stopifnot(
  identical(
    rownames(seuInt[["UMAP3"]]@cell.embeddings), rownames(seuInt[["tSNE3"]]@cell.embeddings)
  )
)

stopifnot(
  identical(
    rownames(seuInt[["UMAP3"]]@cell.embeddings), colnames(spe)
  )
)

reducedDim(spe, "PRECAST_UMAP") <- seuInt[["UMAP3"]]@cell.embeddings
reducedDim(spe, "PRECAST_TSNE") <- seuInt[["tSNE3"]]@cell.embeddings

# Load PRECAST07 annotation
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))



# Plot ----
## UMAP -----
spd_order <- order(spd_anno_df$anno_lab)

### SpD (PRECAST_07) ----
UMAP_spd <- plotReducedDim(
  spe,
  dimred = "PRECAST_UMAP",
  color_by = "PRECAST_07",
  theme_size = 12,
  point_size = 0.1
  # scattermore = TRUE
) +
  scale_color_manual(
    name = "Spatial Domains",
    values = set_names(
      Polychrome::palette36.colors(7)[seq.int(7)],
      unique(spe$PRECAST_07) |> sort()
    )[spd_order],
    labels = spd_anno_df$anno_lab[spd_order]
  ) +
  guides(
    color = guide_legend(title.position = "top", override.aes = list(size = 5))
  ) +
  labs(x = "UMAP1", y = "UMP2") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.justification = "center"
  )

### diagnosis  ----
UMAP_dx <- plotReducedDim(
  spe,
  dimred = "PRECAST_UMAP",
  color_by = "DX",
  theme_size = 12,
  point_size = 0.1
  # scattermore = TRUE
) +
  scale_color_manual(
    name = "Diagnosis",
    values = c(
      "NTC" = "blue",
      "SCZ" = "red"
    )
  ) +
  guides(
    color = guide_legend(title.position = "top", override.aes = list(size = 5))
  ) +
  labs(x = "UMAP1", y = "UMP2") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.justification = "center"
  )

### Paneled plot ----
cowplot::plot_grid(
  UMAP_spd,
  UMAP_dx,
  nrow = 1,
  align = c("h")
)




## TSNE -----
### SpD (PRECAST07) ----
tsne_spd <- plotReducedDim(
  spe,
  dimred = "PRECAST_TSNE",
  color_by = "PRECAST_07",
  theme_size = 12,
  point_size = 0.1
  # scattermore = TRUE
) +
  scale_color_manual(
    name = "Spatial Domains",
    values = set_names(
      Polychrome::palette36.colors(7)[seq.int(7)],
      unique(spe$PRECAST_07) |> sort()
    )[spd_order],
    labels = spd_anno_df$anno_lab[spd_order]
  ) +
  guides(
    color = guide_legend(
      title.position = "top", override.aes = list(size = 5)
    )
  ) +
  labs(x = "t-SNE 1", y = "t-SNE 2") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.justification = "center"
  )

### diagnosis ----
tsne_dx <- plotReducedDim(
  spe,
  dimred = "PRECAST_TSNE",
  color_by = "DX",
  theme_size = 12,
  point_size = 0.1
  # scattermore = TRUE
) +
  scale_color_manual(
    name = "Diagnosis",
    values = c(
      "NTC" = "blue",
      "SCZ" = "red"
    )
  ) +
  guides(
    color = guide_legend(
      title.position = "top",
      override.aes = list(size = 5)
    )
  ) +
  labs(x = "t-SNE 1", y = "t-SNE 2") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.justification = "center"
  )

### Paneled plot ----
tsne_panel <- cowplot::plot_grid(
  tsne_dx,
  tsne_spd,
  nrow = 1,
  align = "h"
)

ggsave(
  here(
    "plots/03_visium_spatial_clustering",
    "tsne_PRECAST_07.pdf"
  ),
  tsne_panel,
  height = 5,
  width = 8
)

# Session info -----
sessioninfo()
