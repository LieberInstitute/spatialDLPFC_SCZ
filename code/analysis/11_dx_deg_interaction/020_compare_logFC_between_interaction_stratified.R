# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

# Load data ----
spd_names <- sprintf("spd%02d", 1:7)

file_stratified <- list.files(
  here("processed-data/rds/09_dx_deg_per_spd"),
  full.names = TRUE
)

file_interaction <- list.files(
  here("processed-data/rds/11_dx_deg_interaction"),
  full.names = TRUE
)

spd_names |>
  set_names() |>
  iwalk(.f = function(.spd, idx) {
    stratified_df <- read_csv(
      grep(.spd, file_stratified,
        ignore.case = TRUE, value = TRUE
      )
    )

    interaction_df <- read_csv(
      grep(.spd, file_interaction,
        ignore.case = TRUE, value = TRUE
      )
    )

    merged_df <- inner_join(
      interaction_df, stratified_df,
      by = c("gene_id" = "ensembl")
    )


    ret_p <- ggplot(merged_df) +
      geom_point(
        aes(
          x = logFC_scz, # Layer stratified analysis
          y = logFC # Interaction model
        )
      ) +
      geom_abline(
        intercept = 0, slope = 1, color = "red",
        linetype = "dotted"
      ) +
      labs(
        x = "LogFC (Layer stratified analysis)",
        y = "LogFC (Interaction model)",
        title = paste0(
          idx,
          " - correlation between stratified and interaction model"
        )
      ) +
      annotate("text",
        x = -Inf, y = Inf,
        label = sprintf(
          "r=%.3f",
          # Correlation of the t-statistics
          cor(
            merged_df$logFC_scz,
            merged_df$logFC,
            use = "pairwise.complete.obs"
          )
        ),
        hjust = -0.1, vjust = 2, size = 5, color = "black"
      )

    ggsave(
      filename = here(
        "plots/11_dx_deg_interaction",
        paste0(
          "scatter_logFC_interaction_stratified_",
          # Remove special characters from the index
          gsub("[^[:alnum:]_]", "_", idx),
          ".pdf"
        )
      ),
      plot = ret_p
    )
  })


# Session info ----
session_info()
