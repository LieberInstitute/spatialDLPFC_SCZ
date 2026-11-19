# Visualize the mouse-to-human ortholog matching steps ----
# A. Gene funnel: tested in CROP-seq -> one-to-one human ortholog -> also
#    tested in each of our rankings
# B. Why genes are dropped: ortholog status by mouse gene-name class
# C. Share of each DEG set kept at each step, per cell type and direction
# D. DEG set size after matching, against fgsea's minSize
# E. Expression of DEGs that are kept vs dropped at the ortholog step
#
# Panels C and D use the 13,977-gene universe shared by the domain-adjusted,
# per-domain and PRS rankings.

# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(sessioninfo)
})


# Config ----
min_size <- 15 # must match 02-fgsea_lee_sets_vs_our_rankings.R
ref_ranking <- "Domain-adjusted"

fld_in <- here("processed-data/21_perturb_seq_enrichment/lee_et_al")
fld_plot <- here("plots/21_perturb_seq_enrichment/lee_et_al")
dir.create(fld_plot, showWarnings = FALSE, recursive = TRUE)

celltype_order <- c(
  "Allcells", "Neuron",
  "Glutamatergic", "GABAergic", "Astro", "Oligo", "Microglia", "Vascular"
)

# Reference categorical palette, fixed slot order
col_blue <- "#2a78d6"
col_orange <- "#eb6834"
col_aqua <- "#1baf7a"
col_yellow <- "#eda100"
col_muted <- "#a3a29b"
ink_secondary <- "#5f5e57"

status_levels <- c("one_to_one", "one_to_many", "many_to_one", "no_ortholog")
status_labels <- c(
  one_to_one = "One-to-one (kept)",
  one_to_many = "One mouse to many human",
  many_to_one = "Many mouse to one human",
  no_ortholog = "No human ortholog"
)
status_cols <- c(
  one_to_one = col_blue,
  one_to_many = col_orange,
  many_to_one = col_aqua,
  no_ortholog = col_yellow
)

theme_set(
  theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#e6e5e0", linewidth = 0.3),
      axis.text = element_text(color = ink_secondary),
      plot.title = element_text(face = "bold", size = 10),
      plot.title.position = "plot",
      legend.position = "top",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold")
    )
)


# Load data ----
ortho_status <- read_csv(
  file.path(fld_in, "lee_ortholog_status.csv"),
  show_col_types = FALSE
)
ortho_121 <- read_csv(
  file.path(fld_in, "lee_ortholog_1to1.csv"),
  show_col_types = FALSE
)
deg_map <- read_csv(
  file.path(fld_in, "lee_deg_mapping.csv"),
  show_col_types = FALSE
)
ours <- readRDS(file.path(fld_in, "our_rankings.rds"))
fgsea_df <- readRDS(file.path(fld_in, "fgsea_lee_sets_vs_our_rankings.rds"))


# A. Gene funnel ----
# Universes of our rankings: the 13,977-gene rankings share one universe
universe_groups <- list(
  "Domain-adjusted, per-domain, PRS" = names(ours$rankings[[ref_ranking]]),
  "NeuN" = names(ours$rankings$NeuN),
  "Neuropil" = names(ours$rankings$Neuropil),
  "PNN" = names(ours$rankings$PNN),
  "Vasc" = names(ours$rankings$Vasc)
)
stopifnot(
  map_lgl(
    ours$meta$ranking[ours$meta$family %in% c("per_domain", "prs")],
    ~ setequal(names(ours$rankings[[.x]]), universe_groups[[1]])
  )
)

n_tested <- nrow(ortho_status)
funnel_df <- bind_rows(
  tibble(step = "Tested in CROP-seq", n = n_tested, stage = "mouse"),
  tibble(step = "One-to-one human ortholog", n = nrow(ortho_121), stage = "mapped"),
  imap(universe_groups, ~ tibble(
    step = paste0("Also tested in ours: ", .y),
    n = sum(ortho_121$human_ensembl %in% .x),
    stage = "ours"
  )) |> list_rbind()
) |>
  mutate(
    step = fct_rev(fct_inorder(step)),
    label = sprintf("%s (%s)", comma(n), percent(n / n_tested, accuracy = 1))
  )

p_a <- ggplot(funnel_df, aes(x = n, y = step, fill = stage)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = label),
    hjust = -0.05, size = 2.8, color = ink_secondary
  ) +
  scale_fill_manual(
    values = c(mouse = col_muted, mapped = col_blue, ours = col_blue),
    guide = "none"
  ) +
  scale_x_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.35))
  ) +
  labs(
    title = "A. Genes kept at each matching step",
    x = "Genes", y = NULL
  ) +
  theme(panel.grid.major.y = element_blank())


# B. Ortholog status by gene-name class ----
status_df <- ortho_status |>
  mutate(
    gene_class = case_when(
      str_detect(mouse_gene, "^Gm[0-9]+$") ~ "Gm (predicted)",
      str_detect(mouse_gene, "Rik[0-9]*$") ~ "Rik (RIKEN cDNA)",
      TRUE ~ "Named genes"
    ),
    gene_class = factor(
      gene_class,
      levels = c("Named genes", "Rik (RIKEN cDNA)", "Gm (predicted)")
    ),
    status = factor(status, levels = rev(status_levels))
  ) |>
  count(gene_class, status)

p_b <- ggplot(status_df, aes(x = n, y = fct_rev(gene_class), fill = status)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.4) +
  scale_fill_manual(
    values = status_cols,
    labels = status_labels,
    breaks = status_levels
  ) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "B. Why genes are dropped",
    subtitle = "Most unmatched genes are unnamed Gm and Rik genes",
    x = "Tested mouse genes", y = NULL
  ) +
  guides(fill = guide_legend(nrow = 2)) +
  theme(panel.grid.major.y = element_blank())


# C. Share of each DEG set kept ----
set_df <- fgsea_df |>
  filter(ranking == ref_ranking) |>
  select(celltype, pert, dir, n_set, n_in_ranking, tested)

retain_df <- deg_map |>
  group_by(celltype, pert, dir) |>
  summarise(
    n_deg_mouse = n(),
    n_mapped = sum(!is.na(human_ensembl)),
    .groups = "drop"
  ) |>
  left_join(set_df, by = c("celltype", "pert", "dir")) |>
  mutate(
    n_in_ranking = replace_na(n_in_ranking, 0L),
    `One-to-one human ortholog` = n_mapped / n_deg_mouse,
    `Also tested in ours` = n_in_ranking / n_deg_mouse
  ) |>
  pivot_longer(
    c(`One-to-one human ortholog`, `Also tested in ours`),
    names_to = "step", values_to = "frac"
  ) |>
  mutate(
    step = factor(step, levels = c("One-to-one human ortholog", "Also tested in ours")),
    celltype = factor(celltype, levels = celltype_order),
    dir = factor(dir, levels = c("down", "up"), labels = c("Down DEGs", "Up DEGs"))
  )

retain_med <- retain_df |>
  group_by(celltype, dir, step) |>
  summarise(frac = median(frac), .groups = "drop")

dodge <- position_dodge(width = 0.6)
p_c <- ggplot(retain_df, aes(x = celltype, y = frac, color = step)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1),
    size = 1.1, alpha = 0.55
  ) +
  geom_point(
    data = retain_med,
    position = dodge,
    shape = 95, size = 7
  ) +
  facet_wrap(~dir, nrow = 1) +
  scale_color_manual(values = c(col_blue, col_orange)) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(
    title = "C. Share of each DEG set kept",
    subtitle = paste(
      "Each point is one perturbation; bars are medians.",
      "Up sets lose most genes at the ortholog step"
    ),
    x = NULL, y = "Share of mouse DEGs kept"
  ) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))


# D. Set size after matching, against minSize ----
size_df <- retain_df |>
  filter(step == "Also tested in ours") |>
  mutate(
    status = if_else(n_in_ranking >= min_size, "Tested", "Below minSize (skipped)"),
    status = factor(status, levels = c("Tested", "Below minSize (skipped)"))
  )

p_d <- ggplot(size_df, aes(x = celltype, y = n_in_ranking + 1, color = status)) +
  geom_hline(
    yintercept = min_size + 1,
    linetype = "dashed", color = ink_secondary, linewidth = 0.4
  ) +
  # Label the threshold in the first facet only
  geom_text(
    data = tibble(
      dir = factor("Down DEGs", levels = levels(size_df$dir)),
      label = paste0("minSize = ", min_size)
    ),
    aes(x = Inf, y = (min_size + 1) * 1.35, label = label),
    inherit.aes = FALSE,
    hjust = 1.05, size = 2.6, color = ink_secondary
  ) +
  geom_point(
    position = position_jitter(width = 0.18, height = 0, seed = 1),
    size = 1.3, alpha = 0.8
  ) +
  facet_wrap(~dir, nrow = 1) +
  scale_color_manual(values = c(col_blue, col_muted)) +
  scale_y_log10(
    breaks = c(1, 11, 101, 1001),
    labels = c("0", "10", "100", "1,000")
  ) +
  labs(
    title = "D. DEG set size after matching",
    subtitle = "Sets below minSize are not tested by fgsea",
    x = NULL, y = "Genes in set (log scale)"
  ) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))


# E. Expression of DEGs kept vs dropped ----
expr_df <- deg_map |>
  mutate(
    mapped = if_else(
      !is.na(human_ensembl),
      "Kept (one-to-one)", "Dropped"
    ),
    mapped = factor(mapped, levels = c("Kept (one-to-one)", "Dropped")),
    dir = factor(dir, levels = c("down", "up"), labels = c("Down DEGs", "Up DEGs"))
  )

expr_n <- expr_df |>
  count(dir, mapped) |>
  mutate(label = paste0("n = ", comma(n)))

p_e <- ggplot(expr_df, aes(x = mapped, y = raw_mean_ctrl + 1e-5, fill = mapped)) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA,
    color = ink_secondary, linewidth = 0.3
  ) +
  geom_text(
    data = expr_n,
    aes(x = mapped, y = 2e-6, label = label),
    inherit.aes = FALSE,
    size = 2.5, color = ink_secondary
  ) +
  facet_wrap(~dir, nrow = 1) +
  scale_fill_manual(values = c(col_blue, col_yellow), guide = "none") +
  scale_y_log10(labels = label_log(digits = 1)) +
  labs(
    title = "E. Expression of kept vs dropped DEGs",
    subtitle = "Dropped up DEGs are barely expressed (all cell types, perturbations)",
    x = NULL, y = "Mean UMI per control cell (log scale)"
  )


# Combine and save ----
# wrap_elements() stops patchwork from aligning axes across rows, which
# otherwise pads panel C with panel A's long y labels
p_all <- wrap_elements(full = p_a | p_b) /
  wrap_elements(full = p_c) /
  wrap_elements(full = p_d | p_e) +
  plot_layout(heights = c(1, 1.15, 1.15)) +
  plot_annotation(
    title = "Mouse-to-human ortholog matching of Lee et al. CROP-seq DEGs",
    subtitle = paste0(
      "babelgene, strict one-to-one. DEGs: fdr < 0.05. ",
      "Panels C and D use the ", comma(length(universe_groups[[1]])),
      "-gene universe of the domain-adjusted, per-domain and PRS rankings"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = ink_secondary, size = 9)
    )
  )

ggsave(
  file.path(fld_plot, "ortholog_matching_summary.pdf"),
  p_all,
  width = 12, height = 12,
  bg = "white"
)


# Session info ----
sessioninfo::session_info()
