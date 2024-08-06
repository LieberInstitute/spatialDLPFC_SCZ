library(escheR)
library(SpatialExperiment)
library(tidyverse)
library(limma)
library(sessioninfo)
library(here)
library(spatialLIBD)



# Create PNN+ category in PNN spots
spe <- readRDS(
  here::here(
    "processed-data", "rds", "02_visium_qc",
    # TODO: rename
    "test_qc_spe_w_spg_N63.rds"
  )
)


finalized_spd <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)


col_data_df <- colData(spe) |>
  data.frame() |>
  left_join(
    finalized_spd,
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)


# Filtering only PNN+ spots
spe$pnn_pos <- spe$spg_PWFA >= 0.05

## Visualize PNN spots ----

spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))

pdf(
  here(
    "plots/dx_deg_PNN_pos",
    paste0("test_spot_plot_pnn_pos.pdf")
  ),
  height = 4, width = 5
)
for (.sample_id in unique(spe$sample_id)) {
  plot(
    make_escheR(spe[, spe$sample_id == .sample_id]) |>
      add_ground("PRECAST_07") |>
      add_fill("pnn_pos") +
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
  )
}
dev.off()

## Proportion of the PNN+ spots ----
col_dat <- colData(spe) |> data.frame()

col_dat |>
  group_by(sample_id) |>
  summarize(
    pnn_prop = sum(pnn_pos) / n()
  ) |>
  left_join(
    metadata(spe)$dx_df |> select(sample_id, dx),
    by = "sample_id"
  ) |>
  ggplot() +
  geom_boxplot(
    aes(x = dx, y = pnn_prop, color = dx)
  ) +
  geom_jitter(
    aes(x = dx, y = pnn_prop),
    size = 0.6,
    alpha = 0.4
  ) +
  theme_minimal()


## Proportion of the PNN+ spots per layer ----
col_dat |>
  group_by(sample_id, PRECAST_07) |>
  summarize(
    prop = sum(pnn_pos) / n()
  ) |>
  left_join(
    metadata(spe)$dx_df |> select(sample_id, dx),
    by = "sample_id"
  ) |>
  ggplot() +
  # geom_point(aes(x = PRECAST_07, y = concord, color = dx, group = dx)) +
  geom_boxplot(aes(x = PRECAST_07, y = prop, color = dx)) +
  scale_x_discrete(
    limits = spd_anno_df$spd[order(spd_anno_df$anno_lab)],
    labels = setNames(spd_anno_df$anno_lab, spd_anno_df$spd)
  ) +
  theme_minimal()

# Filter spots ----

pnn_spe <- spe[, spe$pnn_pos == TRUE]

## Pseudobulk based on individual
sce_pseudo <-
  registration_pseudobulk(
    pnn_spe,
    var_registration = "pnn_pos",
    var_sample_id = "sample_id",
    covars = c("dx", "age", "sex", "lot_num", "slide_id"),
    min_ncells = 10
    # pseudobulk_rds_file = here(
    #   "processed-data", "rds", "layer_spd",
    #   paste0("test_spe_pseudo_", .var, ".rds")
    # )
  )

dx_mod <-
  registration_model(
    sce_pseudo,
    covars = c("age", "sex"),
    var_registration = "dx"
  )

# dx_block_cor <-
#   registration_block_cor(
#     sce_pseudo,
#     registration_model = dx_mod,
#     var_sample_id = "sample_id"
#   )

dx_res <- registration_stats_enrichment(
  sce_pseudo,
  block_cor = NaN,
  covars = c("age", "sex"),
  var_registration = "dx",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

dx_res |>
  mutate(abs_log_FC = abs(logFC_scz)) |>
  arrange(desc(abs_log_FC)) |>
  write_csv(
    here(
      "processed-data/PB_dx_pnn",
      "test_dx_deg_pnn_pos_per_sample.csv"
    )
  )


dx_res |>
  arrange(dfdr_scz) |>
  head()

ggplot(
  dx_res,
  aes(
    x = logFC_scz, y = -log10(fdr_scz),
    color = fdr_scz <= 0.05
  )
) +
  geom_point(alpha = 0.8) +
  # geom_label_repel(
  #   data = impl_gene_df,
  #   aes(label = gene),
  #   force = 2,
  #   nudge_y = 0.1
  # ) +
  # geom_label_repel(
  #   data = out_gene_df,
  #   aes(label = gene),
  #   force = 2,
  #   nudge_y = -0.1
  # ) +
  # labs(
  #   title = paste0(
  #     "Pseudobulk analysis by Dx - ", .spd,
  #     " ( ", n_sig_gene, " sig genes)"
  #   )
  # ) +
  theme_minimal()


# run basic analysis adjusting for age and sex

# Session Info
session_info()
