# # Load packages ----
# suppressPackageStartupMessages({
#   library(here)
#   library(tidyverse)
#   library(SpatialExperiment)
#   library(Seurat)
#   library(PRECAST)
#   library(sessioninfo)
#   library(readxl)
#   library(escheR)
# })
# Configuration ----
# NOTE: change this part in case it woul dbe interested to make this as function.
# n_marker_gene <- 5

# # Prepare data
# ## Load spe object ----
# raw_spe <- readRDS(
#   here::here(
#     "processed-data/rds/02_visium_qc",
#     "qc_spe_wo_spg_N63.rds"
#   )
# )

## 20 Samples ----
# Read Excel file
# sub_samples <- read_xlsx(
#   here::here("code/xenium_panel_design/Xenium_DonorList.xlsx"),
#   col_names = FALSE
# )[, 1:2] |> unlist()


# sub_samples <- sub_samples[!is.na(sub_samples)]

# spe <- raw_spe[, raw_spe$brnum %in% sub_samples]

# Prep supervised marker genes ----
## Deprecated - Select spatialDLPFC marker genes  -----


# stopifnot(all(gene_df$model_type == "enrichment"))
# stopifnot(nrow(gene_df) == 9 * n_marker_gene)

# cat(
#   "NOTE (boyiguo1): ",
#   gene_df$ensembl |> unique() |> length(),
#   " unique gene markers are selected for spatial clustering. \n"
# )

## Deconvo buddies marker genes ----
# gene_df <- read_csv(
#   here(
#     "processed-data",
#     "test_xenium_design_layer_genes_deconvo_buddies.csv"
#   )
# ) |> mutate(
#   ensembl = gene_ensembl
# )


# PRECAST Workflow -----

run_PRECAST_with_genes <- function(spe, gene_names) {
  ## Create Seurat List for PRECAST -----
  seuList <- unique(spe$sample_id) |>
    # seuList <- unique(spe$sample_id)[1:4] |> # Testing
    set_names(unique(spe$sample_id)) |>
    map(.f = function(id) {
      tmp_spe <- spe[, spe$sample_id == id]

      tmp_spe$row <- tmp_spe$array_row
      tmp_spe$col <- tmp_spe$array_col

      CreateSeuratObject(
        counts = as.matrix(counts(tmp_spe)),
        meta.data = data.frame(colData(tmp_spe)),
        project = "PNN"
      )
    })


  ## Set-up PRECAST -----
  print("Start PRECAST")
  set.seed(20240416)
  PRECASTObj <- CreatePRECASTObject(
    seuList = seuList,
    selectGenesMethod = NULL,
    customGenelist = gene_names,
    # Set to 0s to force not removing anything
    premin.spots = 0,
    premin.features = 0,
    postmin.spots = 0,
    postmin.features = 0,
    rawData.preserve = FALSE,
    verbose = TRUE
  )

  # Save mem
  # rm(spe)
  rm(seuList)
  gc(verbose = FALSE)

  ## Add a model setting in advance for a PRECASTObj object.
  ## verbose =TRUE helps outputing the information in the algorithm.
  PRECASTObj <- AddAdjList(
    PRECASTObj,
    type = "fixed_distance",
    platform = "Visium"
  )

  PRECASTObj <- AddParSetting(
    PRECASTObj,
    Sigma_equal = FALSE,
    coreNum = ifelse(
      Sys.getenv("SLURM_NTASKS") == "", # on local machine
      8, # Boyi's computer.
      as.numeric(Sys.getenv("SLURM_NTASKS"))
    ),
    maxIter = 30,
    verbose = TRUE,
    seed = 20240416
  )

  print("NOTE (boyiguo1): Finish PRECAST set-up")


  ## Fit PRECAST model ----
  # k_min <- 2
  # k_max <- 16
  PRECASTObj <- PRECAST(
    PRECASTObj,
    K = 7,
    q = 15 # Arbitary number
  )

  spot_keys <- map(
    PRECASTObj@seulist,
    ~ .x@meta.data$key
  )

  PRECAST_df <- map(
    .x = PRECASTObj@resList,
    .f = function(res, # RECAST result for each k, as a list
                  spot_keys) {
      # Error Prevention: if sample size matches
      if (
        !all(
          sapply(res$cluster, length) == sapply(spot_keys, length)
        )
      ) {
        stop("Unconformatble dimensions")
      }

      ret_vec <- map2(
        .x = res$cluster,
        .y = spot_keys,
        .f = function(vec_cluster, vec_keys) {
          vec_cluster <- sprintf(
            "spd%02d",
            vec_cluster |> as.vector()
          )
          names(vec_cluster) <- vec_keys
          return(vec_cluster)
        }
      ) |> list_c()
      return(ret_vec)
    },
    spot_keys = spot_keys
  )

  ## Organize list to df (n*k) ----
  PRECAST_df_final <- do.call(cbind, PRECAST_df)
  k_clus <- attr(PRECASTObj@resList, "para_settings")$K # Order of K
  colnames(PRECAST_df_final) <- sprintf("PRECAST_%02d_test", k_clus)
  PRECAST_df_final <- PRECAST_df_final |>
    data.frame() |>
    rownames_to_column(var = "key")



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
      PRECAST_df_final,
      by = c("key"),
      relationship = "one-to-one"
    ) |>
    left_join(
      finalized_spd,
      by = c("key"),
      relationship = "one-to-one"
    )

  rownames(col_data_df) <- colnames(spe)
  colData(spe) <- DataFrame(col_data_df)

  return(spe)
}


# result_spe <- run_PRECAST_with_genes(spe, gene_df$ensembl)

# Qualitative Examination ----

# pdf(
#   here(
#     "plots/xenium_panel_design",
#     paste0("test_marker_gene_", n_marker_gene, "_per_layer_deconvo_buddies.pdf")
#   ),
#   height = 4, width = 10
# )
# for (.sample_id in unique(spe$sample_id)) {
#   ggpubr::ggarrange(
#     make_escheR(spe[, spe$sample_id == .sample_id]) |>
#       add_fill("PRECAST_07") +
#       scale_fill_manual(
#         values = set_names(
#           Polychrome::palette36.colors(7)[seq.int(7)],
#           unique(spe[["PRECAST_07"]]) |> sort()
#         ),
#         guide = guide_legend(title = NULL)
#       ) + labs(title = .sample_id),
#     make_escheR(spe[, spe$sample_id == .sample_id]) |>
#       add_fill("PRECAST_07_test") +
#       scale_fill_manual(
#         values = set_names(
#           Polychrome::palette36.colors(7)[seq.int(7)],
#           unique(spe[["PRECAST_07"]]) |> sort()
#         ),
#         guide = guide_legend(title = NULL)
#       ) + labs(title = .sample_id)
#   ) |> print()
# }
# dev.off()




# # Session info ----
# sessioninfo::session_info()
