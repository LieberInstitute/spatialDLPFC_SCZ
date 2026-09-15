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
    set convention Ling et al. used for donor-level scoring; this is the set
    used for our projection onto Visium spots.

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
| `02-project_SNAP_scores.R` | Loads the finalized `spe` and the prepared loadings; for each program, computes a per-spot loading-weighted projection score from log-normalized expression, then z-scores it across all spots ("normalized projection"). Saves a slim `key`-indexed `SNAP_score_df.rds` for merging back into `colData(spe)` in downstream scripts. |
| `03-spot_plot_SNAP_rep_samples.R` | Spatial spot plots (via `escheR`) for the two representative samples, Br8667 (`V13M06-342_D1`, NTC) and Br5973 (`V13M06-343_D1`, SCZ): spot border = spatial domain (`spd_label`), spot fill = normalized SNAP-a / SNAP-n projection score. |
| `04-donor_spd_SNAP_distribution.R` | Donor-level visualizations of the SNAP-a / SNAP-n projection score distribution across donors (boxplot) and across spatial domains (boxplot + donor x domain heatmap), split by diagnosis. |

## Projection method

For a given program (SNAP-a or SNAP-n), let `w_g` be the cNMF loading of gene
`g` from the `top2000` gene set, and `x_gj` be the log-normalized expression
of gene `g` in spot `j` (`logcounts(spe)`, matched by gene symbol). The
per-spot projection score is the loading-weighted average expression:

```
score_j = sum_g(w_g * x_gj) / sum_g(w_g)
```

restricted to the genes present in the Visium panel. This score is then
z-scored across all spots (all 63 donors) to produce the "normalized
projection" used for spatial and donor-level visualization, so that SNAP-a
and SNAP-n scores are on a comparable scale across samples.

## Outputs

* `processed-data/rds/16_integration_ling/SNAP_loadings.rds`
* `processed-data/rds/16_integration_ling/SNAP_score_df.rds`
* `plots/16_integration_ling/spot_plot_SNAPa_rep_samples.pdf`
* `plots/16_integration_ling/spot_plot_SNAPn_rep_samples.pdf`
* `plots/16_integration_ling/SNAP_score_by_donor.pdf`
* `plots/16_integration_ling/SNAP_score_by_spatial_domain.pdf`
* `plots/16_integration_ling/SNAP_score_donor_by_domain_heatmap.pdf`
