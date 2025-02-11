# Load Libray -------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  # library(spatialLIBD)
  library(tidyverse)
  # library(escheR)
  # library(ggpubr)
  library(mclust)
  library(ComplexHeatmap)
  library(sessioninfo)
  
})

# Path --------------------------------------------------------------------
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster")

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform",".rds")
)

# Load data ---------------------------------------------------------------

spe <- readRDS(
  path_PRECAST_int_spe
)

spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )]



# DE comparison -----------------------------------------------------------






col_df <- colData(spe) |> data.frame()

precast_col_names <- colnames(col_df) |> 
 grep("^PRECAST_", x = _, value = TRUE)
# 
# res_list <- list()

order_df <- metadata(spe)$dx_df |> arrange(dx, sample_id)

pdf(
  here("plots/spatial_cluster/assessment/heatmap_spd_prop.pdf"),
  width = 11)
for(.k in precast_col_names){
  # browser()
  col_df[[.k]] <- paste0("spd", col_df[[.k]]) |> factor()
  ret_mat <- table(col_df[[.k]], col_df$sample_id) |> 
    prop.table(margin = 2)
  
  # bottom_ha <- HeatmapAnnotation(dx = order_df$dx)
  
  Heatmap(
    ret_mat[, order_df$sample_id], name = "sample_prop",
    # col = circlize::colorRamp2(c(0,0.5,1),  c("blue", "#EEEEEE", "red")),
          # row_order = levels(col_df[[.k]])#,
          # column_order = order_df$sample_id, 
    # column_order = order_df$sample_id,
    column_split = order_df$dx
    # top_annotation = bottom_ha
    ) |> 
    print()

}

dev.off()
# rm(spe); gc();

# colnames(col_df) |> 
#   grep("^PRECAST_", x = _, value = TRUE) |> 
#   # head() |> 
#   set_names() |>
#   map(.f = ~cat(.x)
#     function(.k){
#           # cat(.k)
#     # browser()
#     # col_df$spd <- factor(paste0("spd", col_df[[.k]]))
#     # ret <- col_df |> 
#     #   group_by(spd) |> 
#     #   group_map(~{
#     #     
#     #     table(.x$sample_id) |> prop.table() |> c() |> t() |> data.frame()
#     #   }) |> 
#     #   list_rbind()
#     # 
#     # rownames(ret) <- levels(col_df$spd)
#     # 
#     # ret
#     
#     # browser()
#     # col_df[,.k] |> head()
#   # }
# )



# Session Info ----
sessioninfo::session_info()