# Load Library --------------------
suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(sessioninfo)
  library(tidyverse)
  library(spatialLIBD)
  library(tidyverse)
  library(scater)
  library(scuttle)
  library(pheatmap)
})


# Load Data -----
## Load NMF Result --------
model <- readRDS(
  here(
    "processed-data/rds/NMF/",
    "test_NMF_all_k100.rds"
  )
)

## Load Spe ---------------
spe <- readRDS(
  here(
    "processed-data/rds/spatial_cluster/PRECAST/",
    # "qc_spe_wo_spg_N63.rds",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)


# spe$dx <- metadata(spe)$dx_df$dx[
#   match(
#     spe$sample_id,
#     metadata(spe)$dx_df$sample_id
#   )
# ]

# spe$sex <- metadata(spe)$dx_df$sex[
#   match(
#     spe$sample_id,
#     metadata(spe)$dx_df$sample_id
#   )
# ]

# spe$age <- metadata(spe)$dx_df$age[
#   match(
#     spe$sample_id,
#     metadata(spe)$dx_df$sample_id
#   )
# ]



# colData(spe) |> colnames()

## Attatch PRECAST to SPE ----
# PRECAST_df <- readRDS(
#   here(
#     "processed-data/rds/spatial_cluster",
#     "PRECAST",
#     "test_clus_label_df_semi_inform_k_2-16.rds"
#   )
# )

# stopifnot(nrow(PRECAST_df) == ncol(spe))

# precast_vars <- grep(
#   "^PRECAST_", colnames(PRECAST_df),
#   value = TRUE
# )
# spe <- spe[, spe$key %in% PRECAST_df$key]
# # raw_spe[, precast_vars] <- PRECAST_df[raw_spe$key, precast_vars]
# col_data_df <- PRECAST_df |>
#   right_join(
#     colData(spe) |> data.frame(),
#     by = c("key"),
#     relationship = "one-to-one"
#   )

# rownames(col_data_df) <- colnames(spe)
# colData(spe) <- DataFrame(col_data_df)

## Attatch NMF result to SPE ----
patterns <- t(model$h) # these are the factors
k <- ncol(patterns)
colnames(patterns) <- paste0("NMF", 1:k)
rownames(patterns) <- spe$key
colData(spe) <- cbind(colData(spe), patterns)

## Create nmf sce -----

colnames_kept <- c(
  "key", "sample_id",
  "PRECAST_07",
  "dx", "sex", "age", "brnum"
)

nmf_sce <- SingleCellExperiment(
  assays = list(nmf = t(patterns)),
  colData = colData(spe)[, colnames_kept]
)


vars <- getVarianceExplained(
  nmf_sce,
  variables = c("PRECAST_07", "dx", "sex", "age"),
  assay.type = "nmf"
)

png(here("plots/NMF/var_explained_all_cell.png"))
plotExplanatoryVariables(vars)
dev.off()

# Deprecated code
# vars_dx_only <- getVarianceExplained(
#   nmf_sce,
#   variables = c("dx"),
#   assay.type = "nmf"
# )
# png(here("plots/NMF/var_explained_all_cell_dx_only.png"))
# plotExplanatoryVariables(vars_dx_only)
# dev.off()


## Pseudo-bulk NMF - sample & spd ----
pb_nmf_sce_spd <- aggregateAcrossCells(
  nmf_sce,
  ids = DataFrame(
    sample_id = nmf_sce$sample_id,
    spd = nmf_sce$PRECAST_07
  ),
  use.assay.type = "nmf",
  statistics = "mean"
)

png(here("plots/NMF/var_explained_spd-PB.png"))
plotExplanatoryVariables(
  pb_nmf_sce_spd,
  variables = c("dx", "sex", "age", "ncells", "spd"),
  assay.type = "nmf"
) |> print()
dev.off()


# Heatmap ----
pb_nmf_mat <- assays(pb_nmf_sce_spd)$nmf |> data.matrix()
tmp_key <- paste(
  pb_nmf_sce_spd$sample_id, 
  pb_nmf_sce_spd$PRECAST_07,
  sep = "_" )



colnames(pb_nmf_mat) <- tmp_key

anna_col_df <- data.frame(
  key = tmp_key,
  brnum = pb_nmf_sce_spd$brnum,
  spd  = pb_nmf_sce_spd$PRECAST_07,
  dx = pb_nmf_sce_spd$dx) |> 
  column_to_rownames("key")

order_index <- order(
  # anna_col_df$spd, 
  anna_col_df$dx,
  anna_col_df$brnum 
)
### Ordered only by dx ----
pb_nmf_mat[, order_index] |>
pheatmap(
  scale = "row",
  cluster_cols = FALSE,
  annotation_col = anna_col_df |> select(dx)
)


### ordered by SPD ----
spd_order_index <- order(
  anna_col_df$spd, 
  anna_col_df$dx,
  anna_col_df$brnum 
)
pb_nmf_mat[, spd_order_index] |>
pheatmap(
  # scale = "row",
  cluster_cols = FALSE,
  annotation_col = anna_col_df |> select(spd, dx)
)


pb_nmf_mat[, spd_order_index] |>
pheatmap(
  # scale = "row",
  cluster_cols = FALSE,
  annotation_col = anna_col_df |> select(spd, dx)
)







## Pseudo-bulk NMF - sample ----
# Deprecated.
# pb_nmf_sce_raw <- aggregateAcrossCells(
#   nmf_sce,
#   ids = nmf_sce$sample_id,
#   use.assay.type = "nmf"
# )
# png(here("plots/NMF/var_explained_sample-PB_raw.png"))
# plotExplanatoryVariables(
#   pb_nmf_sce_raw,
#   variables = c("dx", "sex", "age", "ncells"),
#   assay.type = "nmf"
# ) |> print()
# dev.off()


# pb_nmf_sce <- aggregateAcrossCells(
#   nmf_sce,
#   ids = nmf_sce$sample_id,
#   use.assay.type = "nmf",
#   statistics = "mean"
# )

# png(here("plots/NMF/var_explained_sample-PB.png"))
# plotExplanatoryVariables(
#   pb_nmf_sce,
#   variables = c("dx", "sex", "age", "ncells"),
#   assay.type = "nmf"
# ) |> print()
# dev.off()

p_vec <- apply(assay(pb_nmf_sce, "nmf"),
  MARGIN = 1,
  FUN = function(.nmf) {
    mdl <- lm(.nmf ~ dx + sex + age,
      data = makePerCellDF(pb_nmf_sce)
    )
    summary(mdl)$coef["dxscz", "Pr(>|t|)"]
  }
)

png(here("plots/NMF/sample-PB_qqplot.png"))
qqplot(
  qunif(ppoints(p_vec)),
  p_vec
)
dev.off()



# Spot Plot of NMF Patterns ----
.sample_ordered <- metadata(spe)$dx_df |>
  arrange(dx, sample_id) |>
  pull(sample_id)

for (.fac in paste0("NMF", 1:k)) {
  # .fac <- "NMF1"
  vis_grid_gene(
    spe,
    geneid = .fac,
    pdf_file = here(paste0("plots/NMF/test_spot_plot_", .fac, ".pdf")),
    sample_order = .sample_ordered,
    point_size = 0.8,
    spatial = FALSE
  )
}


# Correlation of Factors ----
mat_NMF_cor <- cor(patterns)
pdf(here("plots/NMF/NMF_corr_mat.pdf"), width = 10, height = 10)
heatmap(mat_NMF_cor, symm = TRUE, scale = "none")
dev.off()

pdf(here("plots/NMF/abs_NMF_corr_mat.pdf"), width = 15, height = 15)
heatmap(abs(mat_NMF_cor), symm = TRUE, scale = "none")
dev.off()





## Scatter plot -----

# t.test(spe$NMF1~spe$dx)
# t.test(spe$NMF2~spe$dx)
pdf(here("plots/NMF/test_nmf_k50_all_spots_density_comp.pdf"))
paste0("NMF", 1:k) |>
  walk(.f = function(nmf_name) {
    (ggplot() +
      geom_density(aes(x = spe[[nmf_name]], group = spe$dx)) +
      labs(title = nmf_name)) |>
      print()
  })
dev.off()

col_df <- colData(spe) |> data.frame()
all_spots_test <- paste0("NMF", 1:k) |>
  map(
    ~ t.test(col_df[[.x]] ~ col_df$dx)
  )

NMF_p <- lapply(all_spots_test,
  FUN = function(.res) {
    # browser()
    .res$p.value
  }
) |>
  unlist() |>
  set_names(paste0("NMF", 1:k))


which(NMF_p == 0)

qqplot(
  -log10(qunif(ppoints(NMF_p))),
  -log10(NMF_p + .Machine$double.xmin)
)

qqplot(tmp, qunif(ppoints(tmp)))

qqnorm(trees$Height)
plot(
  y = sort(trees$Height),
  x = qnorm(ppoints(trees$Height))
)


tmp <- c(
  4.297406246, 0.954763715, 0.763137061, 0.474244034, 0.474244034,
  0.416214305, 0.402158823, 0.371188284, 0.056773325, 0.011020541
)

q_tmp <- qqnorm(tmp)

q_tmp$x^10

q_tmp <- qqplot(tmp, qunif(ppoints(tmp)))

q_tmp$y

-log10(q_tmp$y)

pdf(
  here("plots/NMF/test_nmf_k50_all_spots_registration.pdf"),
  width = 15, height = 15
)

for (spd_var in paste0("PRECAST_", 2:9)) {
  # spd_var <- "PRECAST_9"
  # stopifnot(!all(is.numeric(spe[[spd_var]])))

  dn_mat_spd <- model.matrix(~ factor(spe[[spd_var]]) - 1)

  # kronecker(dn_mat_spd)
  # dn_mat_spd |> str()
  # ggplot() +
  #   geom_density(aes(x = spe$NMF2, group = spe$dx))

  reg_mat <- cor(col_df |> select(starts_with("NMF")), dn_mat_spd)

  heatmap(reg_mat, # Rowv = NA, Colv = NA,
    scale = "none",
    # col = colorRampPalette(c("blue", "white", "red"))(100),
    # symm = TRUE,  # Show the matrix symmetrically
    main = paste0("NMF Registering to ", spd_var),
    xlab = "SpD",
    ylab = "NMF"
  ) |> print()
}




# NMF registration to sample ----------------------------------------------
dn_mat_smp <- model.matrix(~ spe$sample_id - 1)
reg_mat <- cor(col_df |> select(starts_with("NMF")), dn_mat_smp)
heatmap(reg_mat,
  Colv = NA, # Rowv = NA,
  scale = "none",
  # col = colorRampPalette(c("blue", "white", "red"))(100),
  # symm = TRUE,  # Show the matrix symmetrically
  main = "NMF Registration to sample",
  xlab = "Sample",
  ylab = "NMF"
) |> print()

dev.off()






# TODO: Test if the factor is associated with

print(
  "Finish the analysis"
)
# Session Info ------------------------------------------------------------
sessioninfo::session_info()
