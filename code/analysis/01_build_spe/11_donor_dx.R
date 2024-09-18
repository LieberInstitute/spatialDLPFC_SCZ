# Load Pacakges -------------------------------------------------------
library(gtsummary)
library(SpatialExperiment)
library(tidyverse)
library(here)
library(sessioninfo)


# Load SPE Object ---------------------------------------------------------
spe <- readRDS(
  here::here(
    "processed-data",
    "rds", "01_build_spe",
    "raw_spe_wo_SPG_N63.rds"
  )
)


# Organize Data Frame -----------------------------------------------------
demo_df <- metadata(spe)$dx_df

# Error Prevention
if ("Br3942" %in% demo_df$subject) {
  warning("Non-DLPFC Sample *Br3942* is still in the dataset")
  demo_df <- demo_df[which(demo_df$subject != "Br3942"), ]
}
stopifnot(nrow(demo_df) == 63) # Should remove Br3942

demo_df_formated <- demo_df |>
  transmute(
    `Brain Number (ID)` = subject,
    `Date Dxt\'d` = `Dissection Date`,
    Race = race,
    Sex = sex,
    Age = age,
    RIN,
    PMI,
    `PrimaryDx` = factor(
      dx,
      levels = c("ntc", "scz"),
      labels = c("Control", "Schizophrenia")
    )
  )


#  Save Raw Demo Data as CSV ----------------------------------------------
demo_df_formated |>
  write_csv(
    file = here::here(
      "manuscript", "tables",
      "stabl_demo_raw.csv"
    )
  )



# Create Demo Compare Table
demo_df_formated |>
  select(-c(
    `Brain Number (ID)`,
    `Date Dxt\'d`
  )) |>
  tbl_summary(
    by = `PrimaryDx`,
    type = list(
      Age ~ "continuous",
      RIN ~ "continuous",
      PMI ~ "continuous"
    ),
    statistic = list(
      Age ~ "{mean} ({sd})",
      RIN ~ "{mean} ({sd})",
      PMI ~ "{mean} ({sd})"
    )
  ) |>
  # Note: View output on JHPCE
  # as_tibble()
  # Note: save as an xlsx file
  as_hux_xlsx(
    file = here(
      "manuscript", "tables",
      "stabl_demo_compare.xlsx"
    )
  )

# A tibble: 8 × 3
# `**Characteristic**` `**Control**, N = 31` `**Schizophrenia**, N = 32`
# <chr>                <chr>                 <chr>
#   1 Race                 NA                    NA
# 2 CAUC                 31 (100%)             32 (100%)
# 3 Sex                  NA                    NA
# 4 F                    15 (48%)              14 (44%)
# 5 M                    16 (52%)              18 (56%)
# 6 Age                  47 (10)               47 (7)
# 7 RIN                  7.80 (7.35, 8.55)     7.85 (7.18, 8.50)
# 8 PMI                  26 (22, 31)           26 (19, 33)


# Session Info ------------------------------------------------------------
session_info()
