# Load library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(here)
})

# Load data ----
## Load standardized PRS data ----
prs_data <- read_csv(
  here(
    "processed-data/donor_prs",
    "Spatial_DLPFC_SCZ_PRS.csv"
  )
)

## Load diag data ----
diag_data <- metadata(
  readRDS(
    here(
      "processed-data/rds/07_dx_pseudobulk",
      "sce_pseudo_PRECAST07_donor_spd.rds"
    )
  )
)$dx_df

## Merged data ----
ful_dat <- prs_data |>
  right_join(
    diag_data |> select(subject, dx),
    by = c("IID" = "subject")
  )

library(ggbeeswarm)

ggplot(
  ful_dat,
  aes(x = dx, y = PRS, color = dx)
) +
  geom_boxplot(
    outlier.shape = NA, alpha = 0.5, width = 0.3
  ) +
  geom_beeswarm(size = 2, alpha = 0.8) +
  labs(
    x = "Diagnosis Group",
    y = "PRS Score (Normalized)"
    # title = "Swarm Plot of Scaled PRS Scores by Diagnosis Group"
  ) +
  scale_x_discrete(
    labels = c(
      "ntc" = "NTC",
      "scz" = "SCZ"
    )
  ) +
  scale_color_manual(
    values = c(
      "ntc" = "blue",
      "scz" = "red"
    ),
    labels = c(
      "ntc" = "NTC",
      "scz" = "SCZ"
    ),
    guide = "none"
  ) +
  theme_classic(base_size = 6) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 6, face = "bold"),
    axis.title.y = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6)
  )

ggsave(
  here(
    "plots/14_prs_deg",
    "PRS_boxplot_spd.pdf"
  ),
  width = 1.5, height = 2.5, units = "in"
)
