# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(tidyverse)
  library(escheR)
  library(ggpubr)
  library(sessioninfo)
})

# Path ----
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster"
)

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform", ".rds")
)

# Load data -----

spe <- readRDS(
  path_PRECAST_int_spe
)

# Preprocessing ----
## Fetch demo info
spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$sex <- metadata(spe)$dx_df$sex[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$age <- metadata(spe)$dx_df$age[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$lot_num <- sapply(strsplit(spe$sample_id, "_"), function(x) x[1])
spe$slide_id <- sapply(strsplit(spe$lot_num, "-"), function(x) x[2])



# Enrichment test for each SpD ----
vars <- grep("^PRECAST", names(colData(spe)), value = TRUE)

.var <- vars[7]

for (.var in vars) {
  cat(
    "Start building pseudobulk makers for ",
    .var, "\n"
  )

  # Create Pseudobulk data
  spe$spd <- paste0("SpD_", spe[[.var]]) |> factor()
  sce_pseudo <-
    registration_pseudobulk(
      spe,
      var_registration = "spd",
      var_sample_id = "sample_id",
      covars = c("dx", "age", "sex", "lot_num", "slide_id"),
      min_ncells = 10,
      pseudobulk_rds_file = here(
        "processed-data", "rds", "layer_spd",
        paste0("test_spe_pseudo_", .var, ".rds")
      )
    )
}


# # Run PCA
# set.seed(20220423)
# # sce_pseudo <- scater::runMDS(sce_pseudo, ncomponents = 20)
# sce_pseudo <- scater::runPCA(sce_pseudo)





# # if (spe[[.var]] |> unique() |> length() == 2) {
# #   next
# # }



# PNN_modeling_results <- registration_wrapper(
#   sce = spe,
#   var_registration = "spd",
#   var_sample_id = "sample_id",
#   covars = c("sex", "age"),
#   gene_ensembl = "gene_id",
#   gene_name = "gene_name"
# )


## Download spatialDLPFC enrichment result ----

# DLPFC_modeling_results <- fetch_data("spatialDLPFC_Visium_modeling_results")
# DLPFC_layer_anno <- read_csv(
#   here(
#     "code/analysis/visium_spatial_clustering",
#     "bayesSpace_layer_annotations.csv"
#   )
# ) |> filter(
#   grepl("^Sp09", cluster)
# )

# ## Download spatialDLPFC enrichment result ----
# nature_modeling_results <- fetch_data("modeling_results")

# # -------------------------------------------------------------------------


# # .var <- vars[[7]]

# # TODO: write a loop
# for (.var in vars) {
#   print("Start registering ", .var)

#   spe$spd <- paste0("SpD_", spe[[.var]]) |> factor()

#   if (spe[[.var]] |> unique() |> length() == 2) {
#     next
#   }
#   # browser()
#   # TODO: consider to ajdust for age, sex.
#   PNN_modeling_results <- registration_wrapper(
#     sce = spe,
#     var_registration = "spd",
#     var_sample_id = "sample_id",
#     gene_ensembl = "gene_id",
#     gene_name = "gene_name"
#   )

#   # str(PNN_modeling_results)

#   # Extract Enrichment t-statistics ------
#   PNN_t_stats <- PNN_modeling_results$enrichment[
#     ,
#     grep(
#       "^t_stat",
#       colnames(PNN_modeling_results$enrichment)
#     )
#   ]
#   colnames(PNN_t_stats) <- gsub("^t_stat_", "", colnames(PNN_t_stats))




#   # str(DLPFC_modeling_results)

#   DLPFC_t_stats <- DLPFC_modeling_results$enrichment[, grep("^t_stat", colnames(DLPFC_modeling_results$enrichment))]
#   colnames(DLPFC_t_stats) <- gsub("^t_stat_", "", colnames(DLPFC_t_stats))





#   # Registering enrichment results to prev data.
#   cor_layer <- layer_stat_cor(
#     stats = PNN_t_stats,
#     modeling_results = DLPFC_modeling_results,
#     model_type = "enrichment",
#     top_n = 100
#   )

#   colnames(cor_layer) <- DLPFC_layer_anno$layer_combo2[
#     match(colnames(cor_layer), DLPFC_layer_anno$cluster)
#   ]


#   pdf(
#     file = here(
#       paste0(
#         "plots/cluster_anno/test_",
#         "PRECAST_K-", spe[[.var]] |> unique() |> length(),
#         ".pdf"
#       )
#     )
#   )
#   # par(mfrow = c(2,1))
#   layer_stat_cor_plot(
#     cor_layer[, sort(colnames(cor_layer))],
#     max = max(cor_layer)
#   )



#   (ggarrange(
#     spe[, spe$sample_id %in% c("V13M06-342_D1")] |>
#       make_escheR() |>
#       add_fill(var = "spd") +
#       labs(title = "V13M06-342_D1 (NTC)"),
#     spe[, spe$sample_id %in% c("V13M06-343_D1")] |>
#       make_escheR() |>
#       add_fill(var = "spd") +
#       labs(title = "V13M06-343_D1 (SCZ)"),
#     common.legend = TRUE,
#     legend = "bottom"
#   )) |> print()

#   dev.off()
# }


# cat("Job Finished")

# Session Info ----
sessioninfo::session_info()
