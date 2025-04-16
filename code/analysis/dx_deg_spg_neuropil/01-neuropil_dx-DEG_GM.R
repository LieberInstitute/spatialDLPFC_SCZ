# Load library  ----
suppressForeignCheck({
  library(here)
  library(spatialLIBD)
  library(limma)
  library(readxl)
  library(ggrepel)
  library(sessioninfo)
})

# Load data ----
## Load Pseudobulk data ----
raw_pb_neuropil <- readRDS(
  here(
    "processed-data/rds/PB_dx_spg",
    "pseudo_neuropil_pos_donor_spd.rds"
  )
)

## Subset to WM only ----
pb_neuropil <- raw_pb_neuropil[
  ,
  raw_pb_neuropil$registration_variable %in%
    sprintf("spd%02d", c(2, 3, 5, 6, 7))
]

ncol(pb_neuropil)
# [1] 315

# DE analysis ----
dx_mod <-
  registration_model(
    pb_neuropil,
    covars = c("age", "sex"),
    var_registration = "dx"
  )

# 2025-04-14 11:45:46.326476 create model matrix

dx_block_cor <- registration_block_cor(
  pb_neuropil,
  registration_model = dx_mod,
  var_sample_id = "sample_id"
)

# 2025-04-14 11:45:51.214688 run duplicateCorrelation()

dx_res <- registration_stats_enrichment(
  pb_neuropil,
  block_cor = dx_block_cor,
  covars = c("age", "sex", "PRECAST_07"),
  var_registration = "dx",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

# 2025-04-14 11:47:44.999823 computing enrichment statistics
# 2025-04-14 11:47:46.405967 extract and reformat enrichment results

# Exploratory analysis of test restuls ----
## Histogram of p-values ----
hist(
  dx_res$p_value_scz,
  breaks = 100,
  main = "Histogram of p-values",
  xlab = "p-value",
  ylab = "Frequency"
)

## Number of stat sig genes ----

sum(dx_res$fdr_scz < 0.05)
# 145

sum(dx_res$fdr_scz < 0.1)
# 352

## Volcano plot ----
ggplot() +
  geom_point(
    aes(
      x = dx_res$logFC_scz,
      y = -log10(dx_res$p_value_scz),
      color = dx_res$fdr_scz < 0.05
    ),
    alpha = 0.5
  )

dx_res |> arrange(fdr_scz) |> head(10)
dx_res |> arrange(t_stat_scz) |> head(10)
dx_res |> arrange(desc(t_stat_scz)) |> head(10)
dx_res |> arrange(desc(logFC_scz)) |> head(10)


## Save results ----
write_csv(
  dx_res |> arrange(fdr_scz),
  here(
    # TODO: create a results specific folder latter
    "code/analysis/dx_deg_spg_neuropil",
    "neuropil-dx_DEG-GM.csv"
  )
)

write_csv(
  dx_res |> arrange(fdr_scz) |> filter(fdr_scz < 0.05),
  here(
    "code/analysis/dx_deg_spg_neuropil",
    "neuropil-dx_DEG-GM_fdr_05.csv" 
  )
)


## Love lap with 172 dxDEGs ----
dx_deg <- read_csv(
  here(
    "code/analysis/10_dx_deg_adjust_spd",
    "172_prelim_fdr010.csv"
  )
)

# Venn diagram of overlapping genes ----
library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    `all` = dx_deg$ensembl,
    `neuropil+` = dx_res |> filter(fdr_scz < 0.05) |> pull(ensembl)
  ),
  filename = NULL,
  diable_logging = TRUE,
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "Overlap of dx_deg and dx_res Genes"
)

grid.draw(venn.plot)

intersect(
  dx_deg$gene,
  dx_res |> filter(fdr_scz < 0.05) |> pull(gene))



# Annotation
## Find synGO annotated genes ----
### Load SynGO Genes ----
synGO_genes <- read_excel(
  here::here(
    "processed-data/SynGO/syngo_genes.xlsx"
  )
)

### Calcualate odds retios for synGO database enrichment ----

sig_synGO_gene_df <- dx_res |>
  filter(
    fdr_scz < 0.05
  ) |>
  inner_join(
    synGO_genes,
    by = c("ensembl" = "ensembl_id")
  )

sig_synGO_gene_df |> pull(gene)


# Visualization ----
## Fianlized volcano ----
ggplot(
  data = dx_res |> 
  mutate(

  )
  ,
  aes(x = logFC_scz, y = -log10(p_value_scz))
) +
  # Add line for nominal threshold 0.05
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "darkgrey"
  ) +
  # Plot all genes
  geom_point(
    # aes(color = gene_cat),
    size = 2
  ) +
  # Add interested gene labels
  geom_label_repel(
    data = sig_synGO_gene_df, # Add labels last to appear as the top layer
    aes(x = logFC_scz, y = -log10(p_value_scz), label = gene#, color = gene_cat
    ),
    force = 2,
    nudge_y = 0.1,
    size = 3
  ) +
  # Format the legends
  # scale_color_identity(
  #   name = "Significance",
  #   labels = c(
  #     "FDR > 0.10",
  #     "Up-regulated",
  #     "Down-regulated"
  #   )
  # ) +
  labs(
    title = "SCZ-Differentially Expressed Genes among Neuropil+",
    y = "-log10(p-value)",
    x = "log2(Fold Change in SCZ)",
    color = "Significance"
  ) +
  # Format the plot
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top",
    panel.border = element_rect(
      color = "black", fill = NA, size = 1
    )
  )




# Session info ----
session_info()
