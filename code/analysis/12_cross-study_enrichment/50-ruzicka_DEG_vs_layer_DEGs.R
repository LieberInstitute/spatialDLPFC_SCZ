# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggplot2)
  library(sessioninfo)
})

# Load data ----
## Load Ruzicka DEGs -----
ruzicka_deg_list <- readRDS(here(
  "processed-data/rds/12_cross-study_enrichment",
  "ruzicka_cell_type_DEG.rds"
))


ruzicka_deg_list_sig <- ruzicka_deg_list |>
  map(
    ~ .x |>
      filter(ruzicka_sig_gene) |>
      pull(ensembl, name = NULL)
  )

ruzicka_deg_list_up <- ruzicka_deg_list |>
  map(
    ~ .x |>
      filter(
        ruzicka_sig_gene,
        Meta_logFC > 0
      ) |>
      pull(ensembl, name = NULL)
  )

ruzicka_deg_list_down <- ruzicka_deg_list |>
  map(
    ~ .x |>
      filter(
        ruzicka_sig_gene,
        Meta_logFC < 0
      ) |>
      pull(ensembl, name = NULL)
  )

## Load Layer-restricted DEGs -----
layer_restricted_DEG_df <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_restricted_degs_all_spds.csv"
  )
)

overall_deg_list <- layer_restricted_DEG_df |>
  filter(P.Value < 0.05) |>
  group_split(PRECAST_spd) |>
  set_names(
    nm = layer_restricted_DEG_df |>
      distinct(PRECAST_spd) |>
      pull(PRECAST_spd)
  )

up_reg_deg_list <- overall_deg_list |>
  map(~ .x |>
    filter(logFC > 0) |>
    pull(gene_id, name = NULL))

down_reg_deg_list <- overall_deg_list |>
  map(~ .x |>
    filter(logFC < 0) |>
    pull(gene_id, name = NULL))



## Create background gene sets ----
# Intersection of Ruzicka tested genes and layer-tested geness
ruzicka_deg_list |> map(~ .x |> nrow())

layer_restricted_DEG_df |>
  group_by(PRECAST_spd) |>
  summarize(
    n = n()
  )

## Create DEG gene sets -----



# Enrichment analysis ----

up_comp_df <- expand.grid(
  ruzicka_name = names(ruzicka_deg_list_up),
  layer_name = names(up_reg_deg_list),
  stringsAsFactors = FALSE
) |>
  pmap(
    .f = function(ruzicka_name, layer_name) {
      list(
        ruzicka_name = ruzicka_name,
        layer_name = layer_name,
        direction = "up",
        ruzicka_set = list(ruzicka_deg_list_up[[ruzicka_name]]),
        layer_set = list(up_reg_deg_list[[layer_name]]),
        background_set = list(intersect(
          ruzicka_deg_list[[ruzicka_name]] |> pull(ensembl, name = NULL),
          layer_restricted_DEG_df |>
            filter(PRECAST_spd == layer_name) |>
            pull(gene_id)
        ))
      )
    }
  ) |>
  bind_rows()

down_com_df <- expand.grid(
  ruzicka_name = names(ruzicka_deg_list_down),
  layer_name = names(down_reg_deg_list),
  stringsAsFactors = FALSE
) |>
  pmap_dfr(
    .f = function(ruzicka_name, layer_name) {
      list(
        ruzicka_name = ruzicka_name,
        layer_name = layer_name,
        direction = "down",
        ruzicka_set = list(ruzicka_deg_list_down[[ruzicka_name]]),
        layer_set = list(down_reg_deg_list[[layer_name]]),
        background_set = list(intersect(
          ruzicka_deg_list[[ruzicka_name]] |> pull(ensembl, name = NULL),
          layer_restricted_DEG_df |>
            filter(PRECAST_spd == layer_name) |>
            pull(gene_id)
        ))
      )
    }
  ) |>
  bind_rows()





## Fisher's exact test wrapper function ----
enrich_deg_aggreement <- function(
    set_1,
    set_2,
    n_background) {
  # browser()

  # create a list of contigency tables
  # set_1: vector of gene IDs (e.g., DEGs from Ruzicka)
  # set_2: vector of gene IDs (e.g., DEGs from layer-restricted)
  # n_background: total number of genes considered as background

  # Calculate overlap
  n_yes_yes <- length(intersect(set_1, set_2))
  n_yes_no <- length(set_1) - n_yes_yes
  n_no_yes <- length(set_2) - n_yes_yes
  n_no_no <- n_background - (n_yes_yes + n_yes_no + n_no_yes)

  # Build 2x2 contingency table
  cont_tab <- matrix(
    c(n_yes_yes, n_yes_no, n_no_yes, n_no_no),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      "set_1" = c("yes", "no"),
      "set_2" = c("yes", "no")
    )
  )

  fet_res <- fisher.test(cont_tab, alternative = "greater")

  list(
    OR = fet_res$estimate,
    Pval = fet_res$p.value,
    n_overlap = n_yes_yes,
    n_set_1 = length(set_1),
    n_set_2 = length(set_2),
    n_background = n_background
  )
}



# Run wrapper functions ----
up_enrich_long <- up_comp_df |>
  pmap(.f = function(ruzicka_name, layer_name,
                     ruzicka_set, layer_set, background_set,
                     ...) {
    ret_raw <- enrich_deg_aggreement(
      set_1 = ruzicka_set,
      set_2 = layer_set,
      n_background = length(background_set)
    )

    ret <- c(
      ruzicka_name = ruzicka_name,
      layer_name = layer_name,
      direction = "up",
      ret_raw
    )

    return(ret)
  }) |>
  bind_rows()

down_enrich_long <- down_com_df |>
  pmap(.f = function(ruzicka_name, layer_name,
                     ruzicka_set, layer_set, background_set,
                     ...) {
    ret_raw <- enrich_deg_aggreement(
      set_1 = ruzicka_set,
      set_2 = layer_set,
      n_background = length(background_set)
    )

    ret <- c(
      ruzicka_name = ruzicka_name,
      layer_name = layer_name,
      direction = "down",
      ret_raw
    )

    return(ret)
  }) |>
  bind_rows()

# Plot results ----

plot_up_enrich <- up_enrich_long |>
  mutate(
    ruzicka_name = factor(ruzicka_name, levels = names(ruzicka_deg_list_up)),
    layer_name = factor(layer_name, levels = c(
      "SpD07-L1",
      "SpD06-L2/3",
      "SpD02-L3/4",
      "SpD05-L5",
      "SpD03-L6",
      "SpD01-WMtz",
      "SpD04-WM"
    ) |> rev())
  ) |>
  ggplot(aes(
    x = ruzicka_name,
    y = layer_name,
    color = OR,
    size = -log10(Pval)
  )) +
  geom_point() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 1, limits = c(0.5, max(up_enrich_long$OR, na.rm = TRUE)),
    na.value = "lightgray"
  ) +
  labs(
    title = "Ruzicka DEGs vs Layer-restricted DEGs (Up-regulated)",
    x = "Ruzicka DEGs",
    y = "Layer-restricted DEGs",
    color = "Odds Ratio",
    size = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plot_down_enrich <- down_enrich_long |>
  mutate(
    ruzicka_name = factor(ruzicka_name, levels = names(ruzicka_deg_list_down)),
    layer_name = factor(layer_name, levels = c(
      "SpD07-L1",
      "SpD06-L2/3",
      "SpD02-L3/4",
      "SpD05-L5",
      "SpD03-L6",
      "SpD01-WMtz",
      "SpD04-WM"
    ) |> rev())
  ) |>
  ggplot(aes(
    x = ruzicka_name,
    y = layer_name,
    color = OR,
    size = -log10(Pval)
  )) +
  geom_point() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 1, limits = c(0.5, max(down_enrich_long$OR, na.rm = TRUE)),
    na.value = "lightgray"
  ) +
  labs(
    title = "Ruzicka DEGs vs Layer-restricted DEGs (Down-regulated)",
    x = "Ruzicka DEGs",
    y = "Layer-restricted DEGs",
    color = "Odds Ratio",
    size = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Save plots ----
ggsave(
  filename = here("plots/12_cross_study_enrichment",
                  "ruzicka_DEG_vs_layer_DEGs_up.pdf"),
  plot = plot_up_enrich,
  width = 8, height = 3
)
ggsave(
  filename = here("plots/12_cross_study_enrichment",
                  "ruzicka_DEG_vs_layer_DEGs_down.pdf"),
  plot = plot_down_enrich,
  width = 8, height = 3
)


# Session Info ----
session_info()
