# Load libraries ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
})


# Load data ----
## PRS data ----
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

# Caclualte PRS statistics ----
## NTC group ----
ntc_dat <- ful_dat |> filter(dx == "ntc")
prs_median <- median(ntc_dat$PRS, na.rm = TRUE)
prs_iqr <- IQR(ntc_dat$PRS, na.rm = TRUE)

## scz group ----
scz_dat <- ful_dat |> filter(dx == "scz")
scz_prs_median <- median(scz_dat$PRS, na.rm = TRUE)
scz_prs_iqr <- IQR(scz_dat$PRS, na.rm = TRUE)


# Nonparametric test comparing PRS between ntc and scz groups ----
prs_wilcox <- wilcox.test(
  PRS ~ dx,
  data = ful_dat |> filter(dx %in% c("ntc", "scz"))
)
print(prs_wilcox)

# Parametric test (t-test) comparing PRS between ntc and scz groups ----
prs_ttest <- t.test(
  PRS ~ dx,
  data = ful_dat |> filter(dx %in% c("ntc", "scz"))
)
print(prs_ttest)