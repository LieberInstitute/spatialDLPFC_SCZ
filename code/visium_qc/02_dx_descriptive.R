# Descriptive Table -------------------------------------------------------
# TODO: move to an analysis
library(gtsummary)
library(SpatialExperiment)

# expr_meta <- read_csv()
# list.files(here::here("processed-data", "rds", 
#                       "spe", "01_build_spe"))

# TODO: readin spe object
spe <- readRDS(here::here("processed-data", "rds", 
                          "spe", "01_build_spe", "test_spe_raw_36.rds"))

demo_df <- metadata(spe)$dx_df

# Error Prevention
if("Br3942" %in% demo_df$brain_num){
  warning("Br3942 is still in the dataset")
  demo_df <- demo_df[which(demo_df$brain_num != "Br3942"), ]
}
stopifnot(nrow(demo_df) == 63) # Should remove Br3942

# TODO: change variable names

# TODO: save the demo_df to an excel file



# Create Demo Table

demo_df |> 
  # TODO: add more demo vars
  select(age, sex, dx) |>
  tbl_summary(by = dx,
              type = list(age ~ "continuous"),
              statistic = list(age ~ "{mean} ({sd})")
              ) |> 
  gtsummary::as_tibble() # |>
  # gt::gtsave(here("code", "visium_qc", "demo_tab.pdf"))

# TODO: save Supplementary Table 2


# `**Characteristic**` `**ntc**, N = 17` `**scz**, N = 18`
# <chr>                <chr>             <chr>            
# 1 age                  46 (9)            47 (7)           
# 2 sex                  NA                NA               
# 3 F                    9 (53%)           8 (44%)          
# 4 M                    8 (47%)           10 (56%)    
