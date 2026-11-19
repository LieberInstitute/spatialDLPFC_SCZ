# Build the Lee et al. (2026) perturbation gene-set library ----
# - Maps mouse genes to human Ensembl IDs (strict one-to-one orthologs).
# - Defines up and down DEG sets per cell-type grouping x perturbation
#   (fdr < 0.05, direction from the sign of glm_coef).
# - Keeps the tested-gene universe per cell-type grouping x perturbation.
# - Saves a slim human-mapped GLM table for the reverse-direction analysis.

# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(babelgene)
  library(sessioninfo)
})


# Config ----
fdr_cut <- 0.05

fld_out <- here("processed-data/21_perturb_seq_enrichment/lee_et_al")
dir.create(fld_out, showWarnings = FALSE, recursive = TRUE)


# Load data ----
## Full GLM results from the Won lab ----
load(
  here(
    "processed-data/perturb-seq-degs/2_won",
    "R01_12_perts_full_glm_result.rda"
  )
)

glm_df <- full_glm_result |>
  as_tibble() |>
  select(
    celltype, pert = test_gRNAtarget, mouse_gene = testgene,
    n_ctrlcells, n_pertcells,
    glm_coef, glm_p, fdr, raw_mean_ctrl, raw_mean_pert
  )
rm(full_glm_result)
gc()

stopifnot(
  !anyDuplicated(glm_df[, c("celltype", "pert", "mouse_gene")])
)


# Mouse to human orthologs ----
# NOTE: `top = FALSE` keeps every ortholog call so that one-to-many and
# many-to-one mappings can be detected and dropped. With the default
# `top = TRUE`, only the best-supported match is kept, which would hide them.
ortho_raw <- orthologs(
  genes = unique(glm_df$mouse_gene),
  species = "mouse",
  human = FALSE,
  top = FALSE
) |>
  filter(human_ensembl != "")

ortho_121 <- ortho_raw |>
  group_by(symbol) |>
  filter(n_distinct(human_ensembl) == 1) |>
  ungroup() |>
  group_by(human_ensembl) |>
  filter(n_distinct(symbol) == 1) |>
  ungroup() |>
  # A few Ensembl IDs carry more than one human symbol (e.g., readthrough
  # loci). Keep one symbol per pair, preferring the upper-case mouse symbol.
  group_by(symbol, human_ensembl) |>
  arrange(desc(human_symbol == toupper(symbol)), human_symbol, .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  distinct(
    mouse_gene = symbol,
    human_ensembl,
    human_symbol
  )

stopifnot(
  !anyDuplicated(ortho_121$mouse_gene),
  !anyDuplicated(ortho_121$human_ensembl)
)

message(
  "Mouse genes tested: ", n_distinct(glm_df$mouse_gene),
  "; with a one-to-one human ortholog: ", nrow(ortho_121)
)

write_csv(ortho_121, file.path(fld_out, "lee_ortholog_1to1.csv"))

## Ortholog status of every tested mouse gene (for 06-plot_ortholog_matching.R) ----
# - one_to_one: kept
# - one_to_many: the mouse gene maps to more than one human gene
# - many_to_one: the mouse gene maps to one human gene that other mouse genes
#   also map to
# - no_ortholog: no human ortholog with an Ensembl ID in babelgene
ortho_status <- tibble(mouse_gene = unique(glm_df$mouse_gene)) |>
  left_join(
    ortho_raw |>
      group_by(symbol) |>
      summarise(n_human = n_distinct(human_ensembl), .groups = "drop"),
    by = c("mouse_gene" = "symbol")
  ) |>
  mutate(
    status = case_when(
      mouse_gene %in% ortho_121$mouse_gene ~ "one_to_one",
      is.na(n_human) ~ "no_ortholog",
      n_human > 1 ~ "one_to_many",
      TRUE ~ "many_to_one"
    )
  ) |>
  select(mouse_gene, status)

stopifnot(sum(ortho_status$status == "one_to_one") == nrow(ortho_121))
print(count(ortho_status, status))

write_csv(ortho_status, file.path(fld_out, "lee_ortholog_status.csv"))


# Human-mapped GLM table ----
glm_human <- glm_df |>
  inner_join(ortho_121, by = "mouse_gene")

saveRDS(glm_human, file.path(fld_out, "lee_glm_human.rds"))


# Gene sets ----
## Set ID format: <celltype>|<pert>|<dir> ----
deg_df <- glm_df |>
  filter(fdr < fdr_cut) |>
  mutate(dir = if_else(glm_coef > 0, "up", "down")) |>
  left_join(ortho_121, by = "mouse_gene")

# Every DEG with its mapping outcome (for 06-plot_ortholog_matching.R)
write_csv(
  deg_df |>
    select(celltype, pert, dir, mouse_gene, human_ensembl, raw_mean_ctrl),
  file.path(fld_out, "lee_deg_mapping.csv")
)

gene_set_df <- deg_df |>
  filter(!is.na(human_ensembl)) |>
  mutate(set_id = paste(celltype, pert, dir, sep = "|"))

gene_sets <- split(gene_set_df$human_ensembl, gene_set_df$set_id)

## Tested-gene universe per cell-type grouping x perturbation ----
universe_df <- glm_human |>
  mutate(universe_id = paste(celltype, pert, sep = "|"))

universes <- split(universe_df$human_ensembl, universe_df$universe_id)

saveRDS(
  list(
    gene_sets = gene_sets,
    universes = universes,
    fdr_cut = fdr_cut
  ),
  file.path(fld_out, "lee_gene_sets.rds")
)


# Set size summary ----
## Size before and after ortholog mapping ----
set_size_df <- deg_df |>
  group_by(celltype, pert, dir) |>
  summarise(
    n_deg_mouse = n(),
    n_deg_human_1to1 = sum(!is.na(human_ensembl)),
    median_ctrl_mean_mouse = median(raw_mean_ctrl),
    median_ctrl_mean_human_1to1 = median(raw_mean_ctrl[!is.na(human_ensembl)]),
    .groups = "drop"
  ) |>
  left_join(
    glm_df |>
      group_by(celltype, pert) |>
      summarise(
        n_tested_mouse = n(),
        n_pertcells = first(n_pertcells),
        .groups = "drop"
      ),
    by = c("celltype", "pert")
  ) |>
  left_join(
    glm_human |>
      count(celltype, pert, name = "n_tested_human_1to1"),
    by = c("celltype", "pert")
  ) |>
  mutate(frac_deg_mapped = n_deg_human_1to1 / n_deg_mouse) |>
  arrange(celltype, pert, dir)

write_csv(set_size_df, file.path(fld_out, "lee_gene_set_size.csv"))

print(
  set_size_df |>
    group_by(celltype, dir) |>
    summarise(
      n_sets = n(),
      median_n_deg_human = median(n_deg_human_1to1),
      median_frac_mapped = median(frac_deg_mapped),
      .groups = "drop"
    ),
  n = Inf
)


# Sanity check: each knocked-out gene is down in its own perturbation ----
target_check <- glm_df |>
  filter(celltype == "Allcells", mouse_gene == pert) |>
  select(pert, glm_coef, fdr)
print(target_check)


# Session info ----
sessioninfo::session_info()
