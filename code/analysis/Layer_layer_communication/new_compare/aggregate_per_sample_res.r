# Load Packages ----
library(CellChat)
library(here)
library(tidyverse)
library(sessioninfo)

# Load Data -----
rds_files <- list.files(
  here(
    "processed-data/layer_layer_comm/per_sample"
  )
)

tmp <- rds_files |>
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

# Convert list of named vectors to a dataframe
# TODO: debug this part

df <- enframe(tmp) %>%
  as.data.frame()

# Rename columns to match the names of the list elements
colnames(df) <- names(tmp)

# Calcualte the pathway for each object


# Create Union of the pathway
unique(c(names(tmp)))

# Create a matrix to compiling the pathways across all samples contianing all samples


# Viz

# Doing t-test between the group
