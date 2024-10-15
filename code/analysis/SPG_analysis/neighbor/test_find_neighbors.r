# Load packages----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
  library(spatialLIBD)
  library(BayesSpace)
  library(escheR)
})


# Load test data -----
# NOTE: one spe object that has spot calling for SPG channels
sub_spe <- readRDS(
  here(
    "processed-data/rds/PB_dx_spg",
    "test_small_spe_for_neighbor.rds"
  )
)


# This is something that Boyi tested for his other projects
sub_spe$array_row |>
  unique() |>
  length()
sub_spe$array_col |>
  unique() |>
  length()
sub_spe$array_row |> table()
sub_spe$array_col |> table()

neighbors_list <-
  BayesSpace:::.find_neighbors(sub_spe, platform = "Visium")


### Helpfer function ----
# NOTE: extract neighbors of specific SPG channel
which_neighbors <- function(spe, var, return_keys = TRUE) {
  # browser()
  # TODO: edit this part values is just boeleans
  # i <- which(colData(spe)[[var]] %in% values)
  i <- which(colData(spe)[[var]] == TRUE)

  res <- sort(unique(unlist(neighbors_list[i]))) + 1
  res_keys <- spe$key[res]


  if (return_keys == TRUE) {
    return(res_keys)
  } else {
    return(res)
  }
}

neighbor_pnn_pos <- which_neighbors(sub_spe, "pnn_pos", return_keys = TRUE)


sub_spe$pnn_neighbor <- FALSE
sub_spe$pnn_neighbor[sub_spe$key %in% neighbor_pnn_pos] <- TRUE

library(escheR)
make_escheR(sub_spe) |>
add_fill("pnn_pos") |>
add_symbol("pnn_neighbor", size = 0.5) +
scale_shape_manual(
  breaks = c("FALSE","TRUE"),
  values = c(NA, 2),
  limits = c("TRUE")
)
