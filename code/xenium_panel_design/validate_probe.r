# Load Packages ----
library(tidyverse)
library(here)
library(sessioninfo)


# Load Probe file ----
prob_df <- read_xlsx(
  here("code/xenium_panel_design/Xenium_SCZ_ProbeSelection5_SHK_v1.xlsx"),
  sheet = "SHK_Final_300_v1",
  col_names = FALSE
) |> select(-`...6`)
## Remove useless columns ----
colnames(prob_df) <- c("gene", "cell_type", "deg", "layer", "other")

# Check duplicate gene names
stopifnot(!duplicated(prob_df$gene))
stopifnot(nrow(prob_df) == 300)

# Number of Cell type markers
sum(!is.na(prob_df$cell_type))
# [1] 65

# Number of Layer markers
sum(!is.na(prob_df$layer))
# [1] 86

sum(!is.na(prob_df$deg))
# [1] 169

sum(!is.na(prob_df$other))
# [1] 15



# Session Info -----
sessino_info()