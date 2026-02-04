# Load library ----
library(here)
library(tidyverse)


# Load data ----
shizhong_prs <- read_csv(
    here(
        "processed-data/donor_prs",
        "SCZ_PRS_with_dx_and_genetic_PCs.csv"
    )
) |> mutate(
    PRS_shizhong = PRS,
    norm_PRS_shizhong = scale(PRS_shizhong, center = TRUE, scale = TRUE)
)

nikos_prs <- read_csv(
    here(
        "processed-data/donor_prs/",
        "archive/w_error/Spatial_DLPFC_SCZ_PRS.csv"
    )
)

## Merge prs together ----
merged_dat <- shizhong_prs |>
    inner_join(
        nikos_prs |> select(subject = IID, PRS_Nikos = PRS),
        by = c("subject")
    )

cor(merged_dat$norm_PRS_shizhong, merged_dat$PRS_Nikos)
# [1,] 0.752685

ggplot(merged_dat, aes(x = norm_PRS_shizhong, y = PRS_Nikos, color = dx)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = c("ntc" = "blue", "scz" = "red")) +
    annotate("text", x = 2, y = -2, label = "r = 0.753", 
                     hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold") +
    labs(
        x = "Normalized PRS (Shizhong)",
        y = "PRS (Nikos)",
        color = "Diagnosis",
        title = "Comparison of PRS Versions"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        panel.grid.major = element_line(color = "gray90")
    )
