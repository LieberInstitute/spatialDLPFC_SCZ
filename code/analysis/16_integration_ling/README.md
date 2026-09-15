# 16 - Integration with Ling et al. SNAP gene programs

This directory projects the "Suppressed Neuron/Astrocyte Program" (SNAP) gene
loadings reported in Ling et al. (schizophrenia single-nucleus RNA-seq study)
onto our Visium spatial transcriptomics data (`spe`), and visualizes the
resulting spot-level and donor-level scores.

Two programs are provided in the shared reference file:

* **SNAP-a**: Astrocyte cNMF2 program (Ling et al. Suppl. Table 6)
* **SNAP-n**: L5 IT neuron cNMF6 program (Ling et al. Suppl. Table 7)

(Note: these are referred to as "SNAP-a"/"SNAP-p" in early discussion of this
analysis, but the reference file only contains SNAP-a and SNAP-n; "SNAP-n" is
used throughout this directory.)

## Data Inputs

* `Ling2024_SNAP_ready_to_use.xlsx`: Reference gene loadings shared by the
  user, derived from Ling et al. Only the following sheets are used:
  * `SNAPa_loadings_ranked` / `SNAPn_loadings_ranked`: full, descending-ranked
    cNMF gene loadings for each program.
  * `SNAPa_top2000` / `SNAPn_top2000`: top 2000 genes by loading - the gene
    set convention Ling et al. used for donor-level scoring.

  The file also contains `donor_SNAP_scores` and `donor_metadata` sheets,
  which are Ling et al.'s own published per-donor SNAP-a/SNAP-n scores and
  demographics from their snRNA-seq cohort. These are not used anywhere in
  this analysis (our projection is computed directly from our own Visium
  expression data, not from their scores) and are intentionally not loaded.
* `processed-data/rds/01_build_spe/fnl_spe_kept_spots_only.rds`: finalized
  Visium `spe` object (63 donors, PRECAST-07 spatial domains annotated as
  `spd_label`), produced in
  [01_build_spe](/code/analysis/01_build_spe).

## Workflow

| Script | Description |
| --- | --- |
| `01-prepare_SNAP_loadings.R` | Reads the Ling et al. xlsx loadings and saves a tidy `SNAP_loadings.rds` (full + top-2000 gene loadings for SNAP-a/SNAP-n). |
| `02-project_SNAP_scores.R` | Loads the finalized `spe` and the prepared loadings; for each program, computes loading-weighted per-spot projections from both the top-2000 and full gene sets. Saves a slim `key`-indexed `SNAP_score_df.rds` for merging back into `colData(spe)` in downstream scripts. |
| `03-spot_plot_SNAP_rep_samples.R` | Spatial spot plots (via `escheR`) for the two representative samples, Br8667 (`V13M06-342_D1`, NTC) and Br5973 (`V13M06-343_D1`, SCZ): spot border = spatial domain (`spd_label`), spot fill = SNAP-a / SNAP-n projection score (top-2000 and all-genes), grayscale (white = minimum, black = maximum). |
| `04-donor_spd_SNAP_distribution.R` | Donor-level visualizations of the SNAP-a / SNAP-n projection score distribution across donors (boxplot) and across spatial domains (boxplot + donor x domain heatmap), split by diagnosis; produced separately for the top-2000 and all-genes gene sets. |

## Projection method

For a given program (SNAP-a or SNAP-n) and gene set (top-2000 or full), let
`w_g` be the cNMF loading of gene `g`, and `x_gj` be the log-normalized
expression of gene `g` in spot `j` (`logcounts(spe)`, matched by gene symbol),
restricted to the genes present in the Visium panel. The per-spot projection
score depends on whether the matched loadings are all nonnegative:

* **Top-2000 gene sets** (Ling et al.'s ranked-by-loading convention: all
  loadings are nonnegative) use a loading-weighted *average* expression:

  ```
  score_j = sum_g(w_g * x_gj) / sum_g(w_g)
  ```

  Dividing by the matched loading sum makes this a weighted average, so a
  score does not increase simply because a gene set has more genes or a
  larger total loading.

* **Full gene sets** include negative loadings (genes the program pushes
  down, not just up), so normalizing by `sum(w_g)` is not meaningful (the
  denominator can be small, zero, or have an arbitrary sign). For these, the
  raw loading-weighted sum is used instead, without normalization:

  ```
  score_j = sum_g(w_g * x_gj)
  ```

Each of the four program-by-gene-set scores is saved as-is (no further
scaling/normalization): `SNAPa_top_2000_score` and `SNAPn_top_2000_score`
(weighted-average, top-2000 gene sets), and `SNAPa_all_genes_score` and
`SNAPn_all_genes_score` (raw weighted-sum, full gene sets).

Where a full ranked loading list contains a gene symbol more than once, its
loadings are summed before matching it to the Visium panel; the all-gene score
therefore retains the total contribution of every listed loading.

## Outputs

* `processed-data/rds/16_integration_ling/SNAP_loadings.rds`
* `processed-data/rds/16_integration_ling/SNAP_score_df.rds`
* `plots/16_integration_ling/spot_plot_SNAPa_top_2000_rep_samples.pdf`
* `plots/16_integration_ling/spot_plot_SNAPn_top_2000_rep_samples.pdf`
* `plots/16_integration_ling/spot_plot_SNAPa_all_genes_rep_samples.pdf`
* `plots/16_integration_ling/spot_plot_SNAPn_all_genes_rep_samples.pdf`
* `plots/16_integration_ling/SNAP_score_by_donor_top_2000.pdf`
* `plots/16_integration_ling/SNAP_score_by_donor_all_genes.pdf`
* `plots/16_integration_ling/SNAP_score_by_spatial_domain_top_2000.pdf`
* `plots/16_integration_ling/SNAP_score_by_spatial_domain_all_genes.pdf`
* `plots/16_integration_ling/SNAP_score_donor_by_domain_heatmap_top_2000.pdf`
* `plots/16_integration_ling/SNAP_score_donor_by_domain_heatmap_all_genes.pdf`
