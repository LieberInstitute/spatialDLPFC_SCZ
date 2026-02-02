# Load libraries ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
  library(readxl)
})


# Load data ----
## PRS data ----
prs_data <- read_csv(
  here(
    "processed-data/donor_prs",
    "SCZ_PRS_scores_by_threshold_samples64.csv"
  )
)

summary(prs_data)

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

stopifnot(nrow(ful_dat) == 63)

# Identify PRS statistics ----
prs_cols <- colnames(ful_dat)[startsWith(colnames(ful_dat), "p cutoff")]

prs_cols |>
  set_names() |>
  map_dfr(
    .f = function(prs_col) {
      data.frame(
        p_t_test = t.test(ful_dat[[prs_col]] ~ ful_dat$dx)$p.value,
        p_wilcox = wilcox.test(ful_dat[[prs_col]] ~ ful_dat$dx)$p.value
      )
    },
    .id = "PRS_threshold"
  ) |>
  arrange(p_t_test)

# According to Shizhong Han's suggestion to use PRS with smallest p-values,
# `p cutoff=0.001` should be used

ful_dat$PRS <- ful_dat$`p cutoff=0.001`


# Caclualte PRS statistics ----
## NTC group ----
ntc_dat <- ful_dat |> filter(dx == "ntc")
ntc_prs_median <- median(ntc_dat$PRS, na.rm = TRUE)
# [1] -0.00735863
ntc_prs_iqr <- IQR(ntc_dat$PRS, na.rm = TRUE)


## scz group ----
scz_dat <- ful_dat |> filter(dx == "scz")
scz_prs_median <- median(scz_dat$PRS, na.rm = TRUE)
# [1] -0.007162375
scz_prs_iqr <- IQR(scz_dat$PRS, na.rm = TRUE)
# [1] 0.000283835

# Test difference between NTC and SCZ groups ----
## Nonparametric test comparing PRS between ntc and scz groups ----
prs_wilcox <- wilcox.test(
  PRS ~ dx,
  data = ful_dat |> filter(dx %in% c("ntc", "scz"))
)
print(prs_wilcox)

## Parametric test (t-test) comparing PRS between ntc and scz groups ----
prs_ttest <- t.test(
  PRS ~ dx,
  data = ful_dat |> filter(dx %in% c("ntc", "scz"))
)
print(prs_ttest)


# Save organize PRS results to RDS file ----
## Load Genetic PCs data ----
genetic_pcs <- read_xlsx(
  here(
    "processed-data/donor_prs",
    "gPC.xlsx"
  ),
  sheet = "pheno_PC"
)
### Check matching ----
stopifnot(all(ful_dat$"IID" %in% genetic_pcs$BrNum))

# which brain is not included in the genetic PCs data
setdiff(ful_dat$"IID", genetic_pcs$BrNum)
# [1] "Br8207" "Br8433" "Br8492"
# [4] "Br8514" "Br8537" "Br8667"
# [7] "Br8772"

## Merge genetic PCs to ful_dat ----
ful_dat <- ful_dat |>
  left_join(
    genetic_pcs,
    by = c("IID" = "BrNum")
  )

# Check for missing values
stopifnot(nrow(ful_dat) == 63)
stopifnot(!anyNA(ful_dat))

## Check if there's race differences ----
stopifnot(all(
  sort(ful_dat$IID) == sort(diag_data$subject)
))

# Check if donor race is consistent
stopifnot(
  all(ful_dat[order(ful_dat$IID), "Race"] ==
    diag_data[order(diag_data$subject), "race"])
)

## Save the merged file ----
ful_dat |>
  select(
    subject = IID,
    dx,
    PRS,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5
  ) |>
  write_csv(
    here(
      "processed-data/donor_prs",
      "mock_SCZ_PRS_with_dx_and_genetic_PCs.csv"
    )
  )
