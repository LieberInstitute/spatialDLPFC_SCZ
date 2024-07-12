library(tidyverse)
library(here)

# Plot Donor Comparison
donor_df <- readr::read_csv(
  here(
    "raw-data/experiment_info",
    "study_compare.csv"
  )
) |>
  mutate(
    dx = factor(dx, levels = c("SCZ", "NTC")),
    study = factor(study,
      levels = c(
        "Kwon, Guo et al.",
        "Huuki-Myers et al. (2024)",
        "Maynard, Collado-Torres et al. (2021)"
      )
    )
  )

for (tmp_col in c("donor", "sample", "spot")) {
  pdf(
    here(
      "plots/misc",
      paste0("bar_plot_study_compare_", tmp_col, "_num.pdf")
    ),
    height = 3,
    width = 10
  )
  print(
    ggplot(donor_df) +
      geom_bar(aes(y = study, weight = .data[[tmp_col]], fill = dx)) +
      labs(
        x = paste0("# of ", tmp_col, "s"),
        y = ""
      ) +
      scale_fill_manual(
        values = c(
          "NTC" = "blue",
          "SCZ" = "red"
        )
      ) +
      theme_light(base_size = 20)
  )
  dev.off()
}
