# Load Packages ----
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(here)
  library(tictoc)
  library(SSLB)
  library(pheatmap)
  library(sessioninfo)
})

# Load PB Data ----
.var <- "PRECAST_07"
sce <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    paste0("test_spe_pseudo_", .var, ".rds")
  )
)

## Preparing data input for  ----
y_mat <- logcounts(sce) |>
  as.matrix() |>
  t() # sample-by-gene matrix

# Set seed for initial values
set.seed(1)
tic()
crude_mdl <- SSLB(
  Y = y_mat,
  K_init = 10
)
toc()


saveRDS(crude_mdl, "~/SSL_Bicluster.rds")
crude_mdl <- readRDS("~/SSL_Bicluster.rds")

X_SSLB <- crude_mdl$X
K_SSLB <- crude_mdl$K
B_SSLB <- crude_mdl$B


bic_pal <- colorRampPalette(rev(c("#ca0020", "#f4a582", "#ffffff", "#92c5de", "#0571b0")))(100)

cc <- rainbow(nrow(X_SSLB))

rownames(X_SSLB) <- rownames(y_mat)

dx_cc <- rainbow(length(unique(sce$dx))) |> setNames(unique(sce$dx))

dx_cc <- dx_cc[sce$dx]

pheatmap(X_SSLB,
  # Rowv = NULL, Colv = NA,
  border_color = "transparent",
  cluster_cols = FALSE,
  cellheight = 1, cellwidth = 3,
  scale = "column",
  color = bic_pal,
  show_rownames = TRUE,
  # RowSideColors = dx_cc
)

# Session Info ----
sessioninfo::session_info()
