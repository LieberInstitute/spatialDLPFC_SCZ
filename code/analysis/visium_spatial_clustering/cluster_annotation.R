# Load Libray -------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(tidyverse)
  library(escheR)
  library(ggpubr)
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

col_df <- colData(spe) |> data.frame()

# spe_backup <- spe # for debug
# subset_id <- metadata(spe)$dx_df |> group_by(dx) |>
#   slice_head(n=4) |> ungroup() |>
#   pull(sample_id) |>
#     c("V13M06-343_D1", "V13M06-342_D1")
# 
# spe <- spe_backup[ ,spe_backup$sample_id %in% subset_id]


## Down load spatialDLPFC modeling result ----------------------------------

DLPFC_modeling_results <- fetch_data("spatialDLPFC_Visium_modeling_results")
DLPFC_layer_anno <- read_csv(
  here("code/analysis/visium_spatial_clustering",
       "bayesSpace_layer_annotations.csv")
) |> filter(
  grepl("^Sp09", cluster)
)


# -------------------------------------------------------------------------
vars <- grep("^PRECAST", names(col_df), value = TRUE)

#.var <- vars[[7]]

# TODO: write a loop
for(.var in vars){
  
  print("Start registering ", .var)
  
  spe$spd <- paste0("SpD_", spe[[.var]]) |> factor()
  
  if(spe[[.var]] |> unique() |> length() ==2)
    next
  # browser()  
  # TODO: consider to ajdust for age, sex.
  PNN_modeling_results <- registration_wrapper(
    sce = spe,
    var_registration = "spd",
    var_sample_id = "sample_id",
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
  )
  
  # str(PNN_modeling_results)
  
  # Extract Enrichment t-statistics ------
  PNN_t_stats <- PNN_modeling_results$enrichment[,
                                                 grep("^t_stat", 
                                                      colnames(PNN_modeling_results$enrichment)
                                                 )
  ]
  colnames(PNN_t_stats) <- gsub("^t_stat_", "", colnames(PNN_t_stats))
  
  
  
  
  # str(DLPFC_modeling_results)
  
  DLPFC_t_stats <- DLPFC_modeling_results$enrichment[, grep("^t_stat", colnames(DLPFC_modeling_results$enrichment))]
  colnames(DLPFC_t_stats) <- gsub("^t_stat_", "", colnames(DLPFC_t_stats))
  
  # cor_layer <- layer_stat_cor(
  cor_layer <- layer_stat_cor(
    stats = PNN_t_stats,
    modeling_results = DLPFC_modeling_results,
    model_type = "enrichment",
    top_n = 100
  )
  
  colnames(cor_layer) <- DLPFC_layer_anno$layer_combo2[
    match(colnames(cor_layer), DLPFC_layer_anno$cluster)
  ] 
  
  
  pdf(
    file = here(
      paste0(
        "plots/cluster_anno/test_",
        "PRECAST_K-", spe[[.var]] |> unique() |> length(),
        ".pdf"
      )
    )
  )
  # par(mfrow = c(2,1))
  layer_stat_cor_plot(
    cor_layer[, sort(colnames(cor_layer))],
    max = max(cor_layer)
  )
  
  

  (ggarrange(
    spe[, spe$sample_id %in% c("V13M06-342_D1")] |> 
      make_escheR() |> 
      add_fill(var = "spd") + 
      labs(title = "V13M06-342_D1 (NTC)"),
    spe[, spe$sample_id %in% c("V13M06-343_D1")] |> 
      make_escheR() |> 
      add_fill(var = "spd") + 
      labs(title = "V13M06-343_D1 (SCZ)"),
    common.legend = TRUE,
    legend = "bottom"
  )) |>  print()
  
  dev.off()
  
  
  
  
  
  # ggarrange(reg_plot,
  #           sample_plot,
  #           nrow = 2)
  # 
}


cat("Job Finished")

sessioninfo::session_info()



