# Load library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(here)
  library(ggbeeswarm)
})

# Load data ----
## Load standardized PRS data ----
prs_data <- read_csv(
  here(
    "processed-data/donor_prs",
    # "mock_SCZ_PRS_with_dx_and_genetic_PCs.csv"
    "SCZ_PRS_with_dx_and_genetic_PCs.csv"
  )
)

# ## Load diag data ----
# diag_data <- metadata(
#   readRDS(
#     here(
#       "processed-data/rds/07_dx_pseudobulk",
#       "sce_pseudo_PRECAST07_donor_spd.rds"
#     )
#   )
# )$dx_df

# ## Merged data ----
# ful_dat <- prs_data |>
#   right_join(
#     diag_data |> select(subject, dx),
#     by = c("IID" = "subject")
#   )

# Differential test of PRS by diagnosis group ----
## Nonparametric test ----
prs_wilcox <- wilcox.test(
  PRS ~ dx,
  data = prs_data
)
# data:  PRS by dx
# W = 215, p-value =
# 6.777e-05
# alternative hypothesis: true location shift is not equal to 0

## T-test ----
prs_ttest <- t.test(
  PRS ~ dx,
  data = prs_data
)
print(prs_ttest)
# data:  PRS by dx
# t = -4.1498, df = 60.5,
# p-value = 0.0001057



# Make plot ----
ggplot(
  prs_data,
  aes(
    x = dx,
    y = scale(PRS, center = TRUE, scale = TRUE),
    color = dx
  )
) +
  geom_boxplot(
    outlier.shape = NA, alpha = 0.5, width = 0.3
  ) +
  geom_beeswarm(size = 2, alpha = 0.8) +
  labs(
    x = "Diagnosis Group",
    # y = "PRS Score (Normalized)"
    y = "PRS Score"
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
    "norm_PRS_boxplot_spd.pdf"
  ),
  width = 1.5, height = 2.5, units = "in"
)
