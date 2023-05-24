# Descriptive Table -------------------------------------------------------
# TODO: move to an analysis
library(gtsummary)

# expr_meta <- read_csv()

# TODO: readin spe object
spe <- readRDS(here::here("processed-data", "rds", 
                          "spe", "01_build_spe", "spe_raw.rds"))

demo_df <- metadata(spe)

# Create Demo Table
demo_df |> 
  select(age, sex, dx) |>
  tbl_summary(by = dx,
              type = list(age ~ "continuous"),
              statistic = list(age ~ "{mean} ({sd})")
              ) |> 
  gtsummary::as_tibble() # |>
  # gt::gtsave(here("code", "visium_qc", "demo_tab.pdf"))


# # A tibble: 4 × 3
# `**Characteristic**` `**ntc**, N = 4` `**scz**, N = 4`
# <chr>                <chr>            <chr>           
# 1 age                  39.4 (7.0)       42.0 (3.9)      
# 2 sex                  NA               NA              
# 3 F                    2 (50%)          1 (25%)         
# 4 M                    2 (50%)          3 (75%)  
