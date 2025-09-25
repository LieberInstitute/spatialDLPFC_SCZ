# Load packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(tidyverse)
  library(ggbeeswarm)
  library(escheR)
  library(sessioninfo)
})


# Load spe data ----
# use V12F14-053_A1 as example
sub_spe <- readRDS(
  here(
    "processed-data/rds/PB_dx_spg",
    "test_small_spe_for_neighbor.rds"
  )
)

pvalb_index <- which(rowData(sub_spe)$gene_name == "PVALB")

# Validate the gene is correct by cehcking ensemble ID
pvalb_ensembl <- rowData(sub_spe)$gene_id[pvalb_index]
# "ENSG00000100362"
# Seems to be the right one

# five point summary statistics for the gene in this one sample
# V12F14-053_A1
logcounts(sub_spe)[pvalb_index, ] |> summary()
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.0000  0.0000  0.2064  0.0000  3.3038

hist(logcounts(sub_spe)[pvalb_index, ], breaks = 40)

# It seems there's a clear separation between 0 and non-0 expression. Hence, we will use the threshold `logcounts > 0` to call PValb + spots

# Limitation: it is possible that in some spots, there's pvalb expression, just not being picked up due to any techinical reasons. IN this case, we are ignoring those spots here.


# Descriptive Analysis ----
sub_spe$pvalb_logcnt = logcounts(sub_spe)[pvalb_index, ]
sub_spe$pvalb_pos <- logcounts(sub_spe)[pvalb_index, ] > 0
sum(logcounts(sub_spe)[pvalb_index, ] > 0) / nrow(logcounts(sub_spe))
# 0.02156192
# ~2% spots that have PVALB expression

## Per SPD gene expression distribution ----
library(ggridges)
data.frame(
  spd = sub_spe$PRECAST_07,
  pvalb_logcnt = logcounts(sub_spe)[pvalb_index, ]
) |>
  ggplot() +
  geom_density_ridges(
    aes(x = pvalb_logcnt, y = spd, fill = spd),
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = "|", point_size = 3, point_alpha = 1, alpha = 0.7,
  ) +
  scale_y_discrete(
    limits = c(
      "spd07", "spd06", "spd02",
      "spd05", "spd03", "spd01", "spd04"
    )
  ) +
  theme_minimal()

## Prevalence & PB expression ----
# How many (# of spots of each SPD)/percentage(color) spots in each spd have are pvalb_pos

data.frame(
  spd = sub_spe$PRECAST_07,
  pvalb_logcnt = logcounts(sub_spe)[pvalb_index, ]
) |>
  mutate(pvalb_pos = pvalb_logcnt > 0) |>
  group_by(spd) |>
  summarize(
    spd_size = n(),
    pval_prop = sum(pvalb_pos) / n()
  ) |>
  ungroup() |>
  mutate(
    spg_type = "pvalb"
  ) |>
  ggplot() +
  geom_point(
    aes(x = spg_type, y = spd, size = spd_size, color = pval_prop)
  ) + 
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal() +
  scale_y_discrete(
    limits = c(
      "spd07", "spd06", "spd02",
      "spd05", "spd03", "spd01", "spd04"
    )
  )



## Expression Distribution ----
## Spot plot for location ----
make_escheR(sub_spe) |>
  add_fill("pvalb_logcnt") |>
  add_ground("PRECAST_07") +
  scale_shape_manual(
    breaks = c("FALSE", "TRUE"),
    values = c(NA, 3),
    limits = c("TRUE")
  ) +
  scale_fill_gradient(low = "white", high = "black") +
  scale_color_discrete(limits = c(
      "spd07", "spd06", "spd02",
      "spd05", "spd03", "spd01", "spd04"
    ))


# Sesssion Info ----
sessioninfo::session_info()
