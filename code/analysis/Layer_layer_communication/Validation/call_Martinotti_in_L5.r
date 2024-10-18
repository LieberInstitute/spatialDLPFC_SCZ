# srun --pty --x11 --mem=150G --time=5:00:00 --cpus-per-task=4 bash
# Load packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(tidyverse)
  library(ggbeeswarm)
  library(escheR)
  library(sessioninfo)
})

# Load Data ----
## Load Spe Object that has spatial domain label ----
spe <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

PRECAST_df_final <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

col_data_df <- PRECAST_df_final |>
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

## Subset to NTC -----
sub_spe <- spe[, spe$dx == "ntc"]

# Use one sample example
test_spe <- sub_spe[, sub_spe$sample_id == sub_spe$sample_id[1]]


# Scatter plot showing marker genes

plot(
  logcounts(L5_spe)[which(rowData(L5_spe)$gene_name == "RELN"), ],
  logcounts(L5_spe)[which(rowData(L5_spe)$gene_name == "PANTR1"), ]
)



# Subset to only Layer 5 (spd05)
L5_spe <- test_spe[, test_spe$PRECAST_07 == "spd05"]

hist(logcounts(L5_spe)[which(rowData(L5_spe)$gene_name == "RELN"), ], breaks = 20)

L5_spe$log_counts <- logcounts(L5_spe)[which(rowData(L5_spe)$gene_name == "SST"), ]

L5_spe$martinotti <- logcounts(L5_spe)[which(rowData(L5_spe)$gene_name == "RELN"), ] > 0

make_escheR(L5_spe) |>
  add_fill("martinotti") +
  labs(title = "logcounts(RELN) > 0")


L5_spe$martinotti <- logcounts(L5_spe)[which(rowData(L5_spe)$gene_name == "SST"), ] > 1.5

make_escheR(L5_spe) |>
  add_fill("martinotti") +
  labs(title = "logcounts(SST) > 1.5")


# Finding Ligand marker genes -----
which(rowData(L5_spe)$gene_name == "GAD1")
which(rowData(L5_spe)$gene_name == "GAD2")
which(rowData(L5_spe)$gene_name == "SLC6A1")
which(rowData(L5_spe)$gene_name == "SLC6A8")

tmp_df <- data.frame(
  martinotti = logcounts(L5_spe)[which(rowData(L5_spe)$gene_name == "RELN"), ] > 0,
  GAD1 = logcounts(L5_spe)[which(rowData(L5_spe)$gene_name == "GAD1"), ],
  GAD2 = logcounts(L5_spe)[which(rowData(L5_spe)$gene_name == "GAD2"), ]
)

ggplot(tmp_df) +
  geom_boxplot(aes(x = martinotti, y = GAD1))

ggplot(tmp_df) +
  geom_boxplot(aes(x = martinotti, y = GAD2)) +
  geom_beeswarm(aes(x = martinotti, y = GAD2))


# Session Info ----
sessioninfo::session_info()
