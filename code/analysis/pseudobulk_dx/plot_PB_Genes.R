# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(pheatmap)
  library(sessioninfo)
})


# Load Data ----

spe_pb <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

metadata(spe_pb)$dx_df |>
  select(
    sample_id, dx
  ) |>
  arrange(dx) |>
  # pull(sample_id) |>
  column_to_rownames("sample_id")


gene_name <- "MBP"
gene_index <- which(rowData(spe_pb)$gene_name == gene_name)
gene_exp <- logcounts(spe_pb)[gene_index, ]

ge_mat <- data.frame(
  gene_exp,
  spd = spe_pb$PRECAST_07,
  sample_id = spe_pb$sample_id
) |>
  pivot_wider(names_from = spd, values_from = gene_exp) |>
  column_to_rownames("sample_id")

ge_mat |>
  pheatmap(
    cellwidth = 8,
    cellheight = 8,
    border_color = "transparent",
    # cluster_rows = FALSE,
    # labels_row = metadata(spe_pb)$dx_df |>
    #  select(sample_id, dx) |> arrange(dx) |>
    #   pull(sample_id),
    annotation_row = metadata(spe_pb)$dx_df |>
      select(sample_id, dx) |> arrange(dx) |>
      column_to_rownames(var = "sample_id")
  )

ge_mat[metadata(spe_pb)$dx_df |>
  select(sample_id, dx) |>
  arrange(dx) |>
  pull(sample_id), ] |>
  pheatmap(
    cellwidth = 8,
    cellheight = 8,
    border_color = "transparent",
    cluster_rows = FALSE,
    # labels_row = ,
    annotation_row = metadata(spe_pb)$dx_df |>
      select(sample_id, dx) |> arrange(dx) |>
      column_to_rownames(var = "sample_id")
  )

# Doesn't work, because NA valu value after subsetting each column
# ge_mat [metadata(spe_pb)$dx_df |>
#    select(sample_id, dx) |> arrange(dx) |>
#     pull(sample_id), "spd06"] |>
# pheatmap(
#   cellwidth = 8,
#   cellheight = 8,
#   cluster_cols = FALSE,
#   border_color = "transparent",
#   cluster_rows = FALSE,
#   # labels_row = ,
#   annotation_row = metadata(spe_pb)$dx_df |>
#    select(sample_id, dx) |> arrange(dx) |>
#    column_to_rownames(var = "sample_id")
# )

# Scatter plot doesn't work because

# Spagetthi plot ----
data.frame(
  gene_exp,
  spd = spe_pb$PRECAST_07,
  sample_id = spe_pb$sample_id,
  dx = spe_pb$dx
) |> ggplot() +
  geom_boxplot(aes(x = spd, y = gene_exp, color = dx)) +
  scale_x_discrete(limits = c(
    "spd07", "spd06", "spd02",
    "spd05", "spd03", "spd01", "spd04"
  ))


data.frame(
  gene_exp,
  spd = spe_pb$PRECAST_07,
  sample_id = spe_pb$sample_id,
  dx = spe_pb$dx
) |> ggplot() +
  geom_point(
    aes(x = spd, y = gene_exp, color = dx, group = sample_id),
    alpha = 0.5, size = 0.5
  ) +
  geom_line(
    aes(x = spd, y = gene_exp, color = dx, group = sample_id),
    alpha = 0.5
  ) +
  scale_x_discrete(
    limits = c(
      "spd07", "spd06", "spd02",
      "spd05", "spd03", "spd01", "spd04"
    )
  ) +
  labs(title = gene_name)

# ----
gene_df <- data.frame(
  gene_exp,
  colData(spe_pb) |> data.frame()
) |> mutate(
  spd = PRECAST_07,
  age_scaled = scale(age)
)

gene_df$spd <- gene_df$PRECAST_07

crude_mdl <- lm(gene_exp ~ age_scaled + sex, data = gene_df)
mdl <- lm(gene_exp ~ age + sex + I(age^2) + spd, data = gene_df)

gene_df$residual <- mdl$residuals

gene_df |> ggplot() +
  geom_point(
    aes(x = spd, y = residual, color = dx, group = sample_id),
    alpha = 0.5, size = 0.5
  ) +
  geom_line(
    aes(x = spd, y = residual, color = dx, group = sample_id),
    alpha = 0.5
  ) +
  scale_x_discrete(
    limits = c(
      "spd07", "spd06", "spd02",
      "spd05", "spd03", "spd01", "spd04"
    )
  ) +
  labs(title = paste0(gene_name, "_residual"))


gene_df |> ggplot() +
  geom_boxplot(
    aes(x = spd, y = residual, color = dx),
    alpha = 0.5, size = 0.5
  ) +
  # geom_line(
  #   aes(x = spd, y = residual, color = dx, group = sample_id),
  #   alpha = 0.5
  # ) +
  scale_x_discrete(
    limits = c(
      "spd07", "spd06", "spd02",
      "spd05", "spd03", "spd01", "spd04"
    )
  ) +
  labs(title = paste0(gene_name, "_residual"))


## understand the size of the residual
library(ggbeeswarm)
gene_df |> ggplot() +
  # geom_boxplot(
  #   aes(x = spd, y = residual, ),
  #   alpha = 0.5, size = 0.5
  # ) +
  geom_boxplot(
    aes(x = spd, y = residual, color = dx),
    alpha = 0.5, size = 0.5
   ) +
  # geom_quasirandom(aes(x = spd, y = residual, group =dx, color = dx), method = "tukeyDense") +
  scale_x_discrete(limits = c(
    "spd07", "spd06", "spd02",
    "spd05", "spd03", "spd01", "spd04"
  ))




mdl_dx <- lm(gene_exp ~ sex + age + I(age^2) + spd + dx, data = gene_df)

mdl_resid_int <- lm(residual ~ dx*spd , gene_df)



mdl_resid <- lm(residual ~ dx, gene_df)
gene_df |> ggplot() +
geom_boxplot(aes(x = spd, y = mdl_resid$residuals, color = dx))



# cor(gene_df$residual, gene_df$gene_exp)







ggplot() +
  geom_tile(aes(x = sample_id, y = spd, fill = gene_exp))
