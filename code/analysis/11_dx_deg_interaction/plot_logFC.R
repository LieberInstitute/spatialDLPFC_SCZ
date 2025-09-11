# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(sessioninfo)
})

# Load data ----
## Load spd annotation ----
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(
    anno_lab = factor(
      paste0(gsub("spd", "SpD", spd), "-", label),
      levels = c(
        "SpD07-L1/M",
        "SpD06-L2/3",
        "SpD02-L3/4",
        "SpD05-L5",
        "SpD03-L6",
        "SpD01-WMtz",
        "SpD04-WM"
      )
    )
  )

## Load layer-adjusted log FC files ----
adj_deg_df_raw <- read_csv(here(
  "processed-data/rds/10_dx_deg_adjust_spd",
  "dx-deg_PRECAST07.csv"
))

## Load layer-restricted logFC files ----
spd_files <- list.files(
  "processed-data/rds/11_dx_deg_interaction", ,
  pattern = "layer_specific_logFC_.*\\.csv",
  full.names = TRUE
)

spd_deg_list <-
  spd_files |>
  set_names(
    # Annotate SpD with layer labels
    nm = spd_anno_df[
      match(
        str_extract(
          spd_files,
          "(?<=layer_specific_logFC_).*?(?=\\.csv)"
        ),
        spd_anno_df$spd
      ),
      "anno_lab",
      drop = TRUE
    ]
  ) |>
  map(
    ~ read_csv(.x) |>
      # Find the gene symbols
      left_join(
        adj_deg_df_raw |> select(ensembl, gene),
        by = c("gene_id" = "ensembl")
      )
  )

# Idenityfy layer-specific DEGs to highlight ----
# NOTE: Top 5 largest up-reg and down-reg genes of each SpD
## Create df for top layer-specific DEGs ----
unique_genes_nom <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_uniquely_specific_genes_nom_p.csv"
  )
) |>
  mutate(
    spd = ifelse(spd == "SpD07-L1", "SpD07-L1/M", spd)
  )

unique_genes_list <- unique_genes_nom |>
  group_split(spd) |>
  set_names(unique(unique_genes_nom$spd))

layer_specific_gene_df <- unique_genes_list |>
  imap_dfr(.f = function(.x, .y) {
    # browser()
    spd_deg_list[[.y]] |>
      right_join(
        .x,
        by = c("gene_id" = "ENSEMBL")
      )
  })

# Rank based on logFC
select_layer_specific_gene_df_logFC <- bind_rows(
  layer_specific_gene_df |> slice_max(logFC, n = 5, by = spd),
  layer_specific_gene_df |> slice_min(logFC, n = 5, by = spd)
) |>
  # Order the spatial domains
  mutate(spd = factor(spd,
    levels = c(
      "SpD07-L1/M",
      "SpD06-L2/3",
      "SpD02-L3/4",
      "SpD05-L5",
      "SpD03-L6",
      "SpD01-WMtz",
      "SpD04-WM"
    )
  )) |>
  arrange(spd, desc(logFC))

write_csv(
  select_layer_specific_gene_df_logFC,
  here(
    "code/analysis/11_dx_deg_interaction",
    "layer_specific_DEG_logFC_top_5.csv"
  )
)

# Ranked based on t-stat
select_layer_specific_gene_df_t <- bind_rows(
  layer_specific_gene_df |> slice_max(t, n = 5, by = spd),
  layer_specific_gene_df |> slice_min(t, n = 5, by = spd)
) |>
  # Order the spatial domains
  mutate(spd = factor(spd,
    levels = c(
      "SpD07-L1/M",
      "SpD06-L2/3",
      "SpD02-L3/4",
      "SpD05-L5",
      "SpD03-L6",
      "SpD01-WMtz",
      "SpD04-WM"
    )
  )) |> arrange(spd, desc(t))


# Ensembl id or gene symbol of selected genes
gene_ensembl <- select_layer_specific_gene_df_logFC |> pull(gene_id)
gene_symbol <- select_layer_specific_gene_df_logFC |> pull(gene)


# Plot heatmap ----
## Layer-restricted logFC heatmap ----
## subseting the unique_genes_list to include the only genes
specific_logFC_df <- spd_deg_list |>
  imap_dfr(
    ~ .x |>
      filter(gene_id %in% gene_ensembl) |>
      mutate(spd = .y)
  )


spc_heatmap_color <- specific_logFC_df |>
  select(
    # TODO: make it gene symbol in the future,
    gene, spd, logFC
  ) |>
  pivot_wider(
    names_from = spd,
    values_from = logFC
  ) |>
  column_to_rownames(var = "gene") |>
  as.matrix()

spc_heatmap_size <- specific_logFC_df |>
  mutate(
    sig_cat = case_when(
      P.Value > 0.05 ~ "nonsig",
      P.Value < 0.05 & adj.P.Val >= 0.05 ~ "Nominal p < 0.05",
      P.Value < 0.05 & adj.P.Val < 0.05 ~ "FDR < 0.05",
    ),
    size = case_when(
      sig_cat == "nonsig" ~ 0.5,
      sig_cat == "Nominal p < 0.05" ~ 1,
      sig_cat == "FDR < 0.05" ~ 1.5
    )
  ) |>
  select(
    gene, spd, size
  ) |>
  pivot_wider(
    names_from = spd,
    values_from = size
  ) |>
  column_to_rownames(var = "gene") |>
  as.matrix()


spd_order <- c(
  "SpD07-L1/M",
  "SpD06-L2/3",
  "SpD02-L3/4",
  "SpD05-L5",
  "SpD03-L6",
  "SpD01-WMtz",
  "SpD04-WM"
)

spc_heatmap_color <- spc_heatmap_color[gene_symbol, spd_order] |> t()
spc_heatmap_size <- spc_heatmap_size[gene_symbol, spd_order] |> t()


### Create heatmp (gene-by-spd) ----
col_fun <- colorRamp2(
  c(min(spc_heatmap_color), 0, max(spc_heatmap_color)),
  c("blue", "white", "red")
)

spc_ht <- Heatmap(
  spc_heatmap_color,
  name = "logFC",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  # Keep cell boundaries with black lines but hide the heatmap elements
  rect_gp = gpar(col = "black", lwd = 0.5, fill = NA),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.circle(
      x = x, y = y,
      r = unit(0.65 * spc_heatmap_size[i, j], "mm"),
      gp = gpar(fill = col_fun(spc_heatmap_color[i, j]), col = NA)
    )
  },
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6, rot = 45, just = "right"),
  heatmap_legend_param = list(
    title = "logFC",
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 6)
  ),
  show_heatmap_legend = FALSE
)

## Layer-adjusted heatmap ----
adj_deg_df <- adj_deg_df_raw |>
  # subset to gene_name
  filter(ensembl %in% gene_ensembl) |>
  mutate(
    spd = "SpD Adjusted",
    sig_cat = case_when(
      p_value_scz > 0.05 ~ "nonsig",
      p_value_scz < 0.05 & fdr_scz >= 0.05 ~ "Nominal p < 0.05",
      p_value_scz < 0.05 & fdr_scz < 0.05 ~ "FDR < 0.05",
    ),
    size = case_when(
      sig_cat == "nonsig" ~ 0.5,
      sig_cat == "Nominal p < 0.05" ~ 1,
      sig_cat == "FDR < 0.05" ~ 1.5
    )
  )

# Prepare color and size matrix for heatmap
adj_ht_color <- adj_deg_df |>
  select(
    gene,
    `Layer-adjusted` = logFC_scz
  ) |>
  column_to_rownames(var = "gene") |>
  as.matrix()

adj_ht_size <- adj_deg_df |>
  select(
    gene, size
  ) |>
  column_to_rownames(var = "gene") |>
  as.matrix()

stopifnot(identical(dim(adj_ht_color), dim(adj_ht_size)))

# Adjust orientation
adj_ht_color <- adj_ht_color[gene_symbol, , drop = FALSE] |> t()
adj_ht_size <- adj_ht_size[gene_symbol, , drop = FALSE] |> t()

### Create heatmap (gene-by-1-row) ----
col_fun <- colorRamp2(
  c(min(spc_heatmap_color), 0, max(spc_heatmap_color)),
  c("blue", "white", "red")
)

adj_ht <- Heatmap(
  adj_ht_color,
  name = "logFC",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  rect_gp = gpar(col = "black", lwd = 0.5, fill = NA),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.circle(
      x = x, y = y,
      r = unit(0.65 * adj_ht_size[i, j], "mm"),
      gp = gpar(fill = col_fun(adj_ht_color[i, j]), col = NA)
    )
  },
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6, rot = 45, just = "right"),
  show_heatmap_legend = FALSE
)

# Save the plot ----
pdf(
  here(
    "plots/11_dx_deg_interaction",
    "dot_plot_logFC_layer_adjusted_N_specific.pdf"
  ),
  width = 6.8, height = 1.8
)
draw(adj_ht %v% spc_ht)
# , annotation_legend_list = lgd_list)
dev.off()


## With legends ----


spc_ht_legend <- Heatmap(
  spc_heatmap_color,
  name = "logFC",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  # Keep cell boundaries with black lines but hide the heatmap elements
  rect_gp = gpar(col = "black", lwd = 0.5, fill = NA),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.circle(
      x = x, y = y,
      r = unit(0.65 * spc_heatmap_size[i, j], "mm"),
      gp = gpar(fill = col_fun(spc_heatmap_color[i, j]), col = NA)
    )
  },
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6, rot = 45, just = "right"),
  heatmap_legend_param = list(
    title = "logFC",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  ),
  show_heatmap_legend = TRUE
)

adj_ht_legend <- Heatmap(
  adj_ht_color,
  name = "logFC",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  rect_gp = gpar(col = "black", lwd = 0.5, fill = NA),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.circle(
      x = x, y = y,
      r = unit(0.65 * adj_ht_size[i, j], "mm"),
      gp = gpar(fill = col_fun(adj_ht_color[i, j]), col = NA)
    )
  },
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6, rot = 45, just = "right"),
  show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
    title = "logFC",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  )
)

lgd_list <- list(
  # dot size for p-value
  Legend(
    labels = c("Not sig.", "Nominal p < 0.05", "FDR < 0.05"),
    title = "Significance",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    type = "points",
    pch = 16,
    legend_gp = gpar(fill = "black", fontsize = 6),
    size = unit(1:3, "mm"),
    labels_gp = gpar(fontsize = 6)
  )
)

pdf(
  here(
    "plots/11_dx_deg_interaction",
    "dot_plot_logFC_layer_adjusted_N_specific_legends.pdf"
  ),
  width = 8, height = 1.8
)
draw(adj_ht_legend %v% spc_ht_legend, annotation_legend_list = lgd_list)
dev.off()
