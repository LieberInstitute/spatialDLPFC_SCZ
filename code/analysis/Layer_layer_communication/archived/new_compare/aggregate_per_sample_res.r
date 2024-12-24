# Load Packages ----
library(CellChat)
library(here)
library(tidyverse)
library(sessioninfo)

# Load Data -----
# Load diagnosis from meta data
# Merge dx information into df
tmp_spe <- readRDS(
  here("processed-data/rds/layer_spd", "test_spe_pseudo_PRECAST_07.rds")
)


# Process cellchat per-sample objects ----
rds_files <- list.files(
  here(
    "processed-data/layer_layer_comm/per_sample"
  )
)

stopifnot(length(rds_files) == 63)

pathway_df_long <- rds_files |>
  set_names() |>
  map_dfr(
    .f = function(file_name) {
      cellchat <- readRDS(
        here(
          "processed-data/layer_layer_comm/per_sample",
          file_name
        )
      )
      res <- apply(cellchat@netP$prob, c(3), sum) |>
        enframe() |>
        mutate(sample_id = str_remove(file_name, "\\.rds"))

      return(res)
    }
  )

pathway_df_long <- pathway_df_long |>
  left_join(
    metadata(tmp_spe)$dx_df |> select(sample_id, dx),
    by = "sample_id"
  )

# Convert list of named vectors to a dataframe
pathway_df <- pathway_df_long |>
  pivot_wider(id_cols = "sample_id", names_from = "name", values_from = "value")

pathway_df <- pathway_df |>
  left_join(
    metadata(tmp_spe)$dx_df |> select(sample_id, dx),
    by = "sample_id"
  )

# Viz the pathway groups
pathway_df_long |>
  ggplot() +
  geom_violin(
    aes(x = dx, y = value, color = dx)
  ) +
  geom_jitter(
    aes(x = dx, y = value, color = dx)
  ) +
  facet_wrap("name", scales = "free_y")


# NOTE:  how to do with the NA values
sum(is.na(pathway_df)) / (dim(pathway_df) |> prod())

(is.na(pathway_df) |> colSums()) / nrow(pathway_df)

t.test(Glutamate ~ dx, pathway_df)
wilcox.test(RA ~ dx, pathway_df)


# Viz

# Doing t-test between the group
