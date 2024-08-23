# Load Pacakges -----
library(gtsummary)
library(tidyverse)
library(here)
library(readxl)
library(sessioninfo)


# Load donor ----
sub_samples <- read_xlsx(
  here::here("code/xenium_panel_design/Xenium_DonorList_Edit.xlsx"),
  col_names = FALSE
)[, 1:2] |> unlist()
sub_samples <- sub_samples[!is.na(sub_samples)]
stopifnot(length(sub_samples) == 24)

# Load PB data to find demo informaiton ----
spe_pb <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)



# Organize Data Frame -----
demo_df <- metadata(spe_pb)$dx_df |>
  filter(subject %in% sub_samples)
stopifnot(nrow(demo_df) == 24)

# Error Prevention
if ("Br3942" %in% demo_df$subject) {
  warning("Non-DLPFC Sample *Br3942* is still in the dataset")
  demo_df <- demo_df[which(demo_df$subject != "Br3942"), ]
}

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
  as_tibble()
  # Note: save as an xlsx file
  # as_hux_xlsx(
  #   file = here(
  #     "manuscript", "tables",
  #     "stabl_demo_compare.xlsx"
  #   )
  # )

# A tibble: 8 × 3
#   `**Characteristic**` `**Control**  \nN = 12` `**Schizophrenia**  \nN = 12`
#   <chr>                <chr>                   <chr>                        
# 1 Race                 NA                      NA                           
# 2 CAUC                 12 (100%)               12 (100%)                    
# 3 Sex                  NA                      NA                           
# 4 F                    7 (58%)                 6 (50%)                      
# 5 M                    5 (42%)                 6 (50%)                      
# 6 Age                  47 (10)                 49 (8)                       
# 7 RIN                  7.92 (0.73)             8.04 (0.61)                  
# 8 PMI                  28 (11)                 27 (9)    


# Session Info -----
session_info()
