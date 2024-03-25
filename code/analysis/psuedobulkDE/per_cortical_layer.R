# Load Libray -------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(sessioninfo)
  library(tidyverse)
  library(ggrepel)
})



# Mem requested
# TODO: adjust
# First try 20G


# Config -----------------------------------------------------------
## Path Config -------------------------------------------------------------
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster"
)

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform", ".rds")
)

# Load data ---------------------------------------------------------------
spe <- readRDS(
  path_PRECAST_int_spe
)

# spe_backup <- spe # for debug
# NOTE: For testing only
# subset_id <- metadata(spe)$dx_df |> group_by(dx) |>
#   slice_head(n=4) |> ungroup() |>
#   pull(sample_id)
#
# spe <- spe_backup[ ,spe_backup$sample_id %in% subset_id]

## Add dx informaiton to the dataset
spe$dx <- metadata(spe)$dx_df$dx[
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


spe$sex <- metadata(spe)$dx_df$sex[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]



# Sample Level Bulk DE ----------------------------------------------------
spe_bulk <-
  registration_pseudobulk(
    spe,
    var_registration = "dx",
    var_sample_id = "sample_id",
    min_ncells = 10
  )

covars <- c("age", "sex")

registration_mod <-
  registration_model(spe_bulk,
    covars = covars
  )


block_cor <-
  registration_block_cor(
    spe_bulk,
    registration_model = registration_mod
  )

results_enrichment <-
  registration_stats_enrichment(
    spe_bulk,
    block_cor = block_cor,
    covars = covars,
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
  )

saveRDS(
  results_enrichment,
  file = here("processed-data/rds/pseudo_bulk/sample_bulk_DE.rds")
)

results_enrichment <- readRDS(here("processed-data/rds/pseudo_bulk/sample_bulk_DE.rds"))


# Save top 100 up/down-reg genes
rbind(
  results_enrichment |> arrange(logFC_scz) |> slice_head(n = 100),
  results_enrichment |> arrange(logFC_scz) |> slice_tail(n=100)
) |>
write.csv(
  "~/bulk_DEG_200.csv"
)


sum(results_enrichment$fdr_ntc <= 0.05)
sum(results_enrichment$p_value_ntc <= 0.05)

hist(results_enrichment$p_value_ntc) # Histogram of p-value


# Volcano plot
# TODO: interpret and confirm the directionality

out_gene_df <- results_enrichment |>
  arrange(fdr_scz) |>
  slice_head(n = 1)

impl_gene_df <- results_enrichment |>
  filter(gene %in% c(
    "PVALB",
    "NOS1",
    "SST",
    "CHODL",
    "GRIN2A",
    "SV2A",
    "DLG4",
    "C4A",
    "C3"
  )) |>
  select(ensembl, gene, ends_with("scz"))


ggplot(
  results_enrichment,
  aes(x = logFC_scz, y = -log10(fdr_scz))
) +
  geom_point() +
  geom_label_repel(
    data = impl_gene_df, # Add labels last to appear as the top layer
    aes(label = gene),
    force = 2,
    nudge_y = 0.1
  ) +
  geom_label_repel(
    data = out_gene_df,
    aes(label = gene),
    color = "red",
    force = 2,
    nudge_y = -0.1
  ) +
  labs(
    title = "Bulk analysis pooling all spots"
  ) +
  theme_minimal()

results_enrichment[which.min(results_enrichment$fdr_scz), ]

results_enrichment |>
  filter(gene %in% c(
    "PVALB",
    "NOS1",
    "SST",
    "CHODL"
  )) |>
  select(ensembl, gene, ends_with("scz"))


results_enrichment |>
  arrange(p_value_ntc) |>
  slice_head(n = 10)

results_enrichment |>
  arrange(desc(abs(logFC_scz))) |>
  slice_head(n = 10)


## Implication Gene Plotting ---------
data.frame(
  dx = spe_bulk$dx,
  PVALB_log = logcounts(spe_bulk)["ENSG00000100362", ],
  NOS1_log = logcounts(spe_bulk)["ENSG00000089250", ],
  SST_log = logcounts(spe_bulk)["ENSG00000157005", ],
  CHODL_log = logcounts(spe_bulk)["ENSG00000154645", ]
) |>
  pivot_longer(
    cols = ends_with("_log")
  ) |>
  ggplot() +
  geom_violin(aes(x = dx, y = value)) +
  facet_wrap(vars(name), scales = "free")

# Pseudo-bulk DE ----------------------------------------------------------

## Subset Layers -----------------------------------------------------------
# TODO: examine every PRECAST K
col_names <- colData(spe) |> names()
spd_var <- grep("^PRECAST", col_names, value = TRUE)

for (.spd_var in spd_var) {
  # .spd_var <- spd_var[[1]]



  spe$spd <- paste0("SpD_", spe[[.spd_var]]) |>
    factor()

  spd_list <- spe$spd |>
    levels() |>
    map(
      .f = function(.spd_lvl, spe = spe) {
        spe[, spe$spd == .spd_lvl]
      },
      spe = spe
    )

  # names(spd_list) <-

  ## Create Pseudobulk data --------------
  # TODO: make the following code to iterate over the number of spds
  # spe_sub <- spd_list[[1]]
  pdf(
    here(
      paste0("plots/PB_DE/test_", .spd_var, ".pdf")
    )
  )
  spd_enrich_res <- spd_list |>
    map(
      .f = function(spe_sub) {
        # browser()

        # spe_sub <- spd_list[[3]]
        # THis is taking the sum. The logcounts is CPM scale.
        sce_pseudo <-
          registration_pseudobulk(
            spe_sub,
            var_registration = "dx",
            var_sample_id = "sample_id",
            min_ncells = 10
          )

        covars <- c("age", "sex")
        registration_mod <-
          registration_model(
            sce_pseudo,
            covars = covars
          )

        block_cor <-
          registration_block_cor(
            sce_pseudo,
            registration_model = registration_mod
          )

        results_enrichment <-
          registration_stats_enrichment(
            sce_pseudo,
            block_cor = block_cor,
            covars = covars,
            gene_ensembl = "gene_id",
            gene_name = "gene_name"
          )

        # browser()
        # ggplot(results_enrichment) +
        #   geom_point(
        #     aes(
        #       x = logFC_scz,
        #       y = -log10(fdr_scz))
        #   ) +
        #   title()
        impl_gene_df <- results_enrichment |>
          filter(gene %in% c(
            "PVALB",
            "NOS1",
            "SST",
            "CHODL",
            "GRIN2A",
            "SV2A",
            "DLG4",
            "C4A",
            "C3"
          )) |>
          select(ensembl, gene, ends_with("scz"))

        (ggplot(
          results_enrichment,
          aes(x = logFC_scz, y = -log10(fdr_scz))
        ) +
          geom_point() +
          geom_label_repel(
            data = impl_gene_df, # Add labels last to appear as the top layer
            aes(label = gene),
            force = 2,
            nudge_y = 0.1
          ) +
          labs(
            title = paste0(
              "PB-", .spd_var,
              "-",
              spe_sub$spd |> unique() |> as.character()
            )
          ) +
          theme_minimal()) |>
          print()

        remain_gene_ids <- intersect(
          c(
            "ENSG00000100362",
            "ENSG00000089250",
            "ENSG00000157005",
            "ENSG00000154645"
          ),
          rowData(sce_pseudo)$gene_id
        )

        ge_df <- remain_gene_ids |>
          map_dfc(~ logcounts(sce_pseudo)[.x, ])


        #
        # (data.frame(
        #   dx = sce_pseudo$dx,
        #   PVALB_log = logcounts(sce_pseudo)["ENSG00000100362",],
        #   NOS1_log = logcounts(sce_pseudo)["ENSG00000089250",],
        #   SST_log = logcounts(sce_pseudo)["ENSG00000157005",]#,
        #   # CHODL_log = logcounts(sce_pseudo)["ENSG00000154645",]
        # ) |>

        if (ncol(ge_df) < 1) {
          (
            ggplot() +
              labs(
                title = paste0(
                  "PB-", .spd_var,
                  "-",
                  spe_sub$spd |> unique() |> as.character()
                )
              )
          ) |> print()
        } else {
          colnames(ge_df) <- paste0(
            rowData(sce_pseudo)[remain_gene_ids, "gene_name"],
            "_log"
          )
          ge_df <- cbind(
            ge_df,
            dx = sce_pseudo$dx
          )

          (
            ge_df |>
              pivot_longer(
                cols = ends_with("_log")
              ) |>
              ggplot() +
              geom_violin(aes(x = dx, y = value)) +
              facet_wrap(vars(name), scales = "free") +
              labs(
                title = paste0(
                  "PB-", .spd_var,
                  "-",
                  spe_sub$spd |> unique() |> as.character()
                )
              )
          ) |> print()
        }
        results_enrichment
      }
    )

  dev.off()

  names(spd_enrich_res) <- spe$spd |> levels()

  saveRDS(
    spd_enrich_res,
    here(
      "processed-data/rds/pseudo_bulk",
      paste0("test_", .spd_var, ".rds")
    )
  )
}







# Session info ------------------------------------------------------------
sessioninfo::session_info()
