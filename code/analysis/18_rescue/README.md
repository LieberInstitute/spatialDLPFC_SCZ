# 18 - RESCUE cell type deconvolution

This directory merges spot-level cell type proportions estimated by RESCUE
(a reference-based deconvolution method) into our finalized Visium `spe`.
It saves the result as a new `spe` object, documents how well the merge
lines up, and makes preliminary visualizations of cell type composition in
the representative samples and across spatial domains and donors.

## Data Inputs

* `processed-data/rds/RESCUE_res/spot_level_results.csv`: RESCUE
  spot-level output (one row per spot). Contains spot metadata (`spot_id`,
  `sample_id`, `barcode`, `array_row`, `array_col`, `domain`,
  `pxl_{row,col}_in_fullres`, `nUMI`, `unattributed_ratio`) and 31
  proportion columns: 30 reference cell types (`Astro`, `EndoMural_01-02`,
  `Excit_01-15`, `Inhib_01-06`, `Micro`, `Oligo_01-03`, `OPC`) plus
  `Unattributed Expression`. Together, the 31 columns sum to 1 in every spot.
* `processed-data/rds/01_build_spe/fnl_spe_kept_spots_only.rds`: finalized
  Visium `spe` (63 donors, PRECAST-07 spatial domains in `spd_label`),
  produced in [01_build_spe](/code/analysis/01_build_spe).

## Workflow

| Script | Description |
| --- | --- |
| `01-merge_RESCUE_to_spe.R` | Runs the merge QC checks (below) and stops on any misalignment. Adds the RESCUE columns to `colData(spe)`, matched on `key` (= RESCUE `spot_id`), and saves the new `spe` plus a slim, `key`-indexed data frame. |
| `02-spot_plot_RESCUE_rep_samples.R` | For the representative samples Br8667 (`V13M06-342_D1`, NTC) and Br5973 (`V13M06-343_D1`, SCZ): `escheR` spot plots for every broad class and every fine cell type (fill = proportion, white = 0 to black = max; border = `spd_label`), plus each sample's mean composition per spatial domain (broad-class stacked bars, with and without unattributed; fine cell type heatmap). |
| `03-donor_spd_RESCUE_composition.R` | Donor-level summaries from the slim data frame (no `spe` load needed): broad-class composition by spatial domain (diagnoses pooled), by spatial domain x diagnosis, and per donor stratified by spatial domain (each with and without unattributed), fine cell type x domain heatmaps (absolute and row-scaled), and boxplots of donor x domain mean proportions by diagnosis. |
| `rescue_plot_utils.R` | Shared helpers sourced by `02`/`03`: broad-class mapping, order and colors, the conditional-composition helper, and the representative sample IDs. |

Each `.R` script has a matching `.sh` SLURM wrapper (JHPCE, `conda_R/devel`).

### Conventions

* **New `colData` columns:** each cell type is stored as `RESCUE_<cell type>`
  (`Unattributed Expression` becomes `RESCUE_Unattributed`). RESCUE's own
  `nUMI` and `unattributed_ratio` are kept as `RESCUE_nUMI` and
  `RESCUE_unattributed_ratio`. The 31 proportion column names are listed in
  `metadata(spe)$RESCUE_celltype_cols`.
* **Proportions are kept as RESCUE reports them**, so `Unattributed` is
  included and is *not* renormalized away. Plots show it as its own
  (gray) category.
* **Broad classes** are the sums of fine types (e.g. `Excit` =
  `Excit_01` + ... + `Excit_15`).
* **Compositional (stacked bar) plots come in two versions:**
  * *with unattributed* (no file suffix): proportions as RESCUE reports
    them; `Unattributed` is shown in neutral gray.
  * *excluding unattributed* (`_excl_unattributed` suffix): the
    conditional composition of attributed signal. `Unattributed` is dropped
    and the remaining classes are renormalized to sum to 1 within each
    group (donor x domain for donor-level plots, sample x domain for the
    representative samples). The renormalization is applied to the
    *averaged* proportions, i.e. `mean(p_ct) / mean(1 - p_unattributed)`.
    This avoids inflating spots with almost no attributed signal, as
    per-spot renormalization would. Donor-level conditional compositions are
    then averaged across donors.
* **Broad-class colors:** Astro `#4DAF6E`, EndoMural `#E7609E`, Micro
  `#5E3C99`, Oligo `#E07B39`, OPC `#D4B243`, Excit `#2E6DB4`, Inhib
  `#E05A42`, Unattributed `#b5b4ae` (gray). Inhib and Oligo are close in
  hue (OKLab ΔE ≈ 7), so the stacking order (Excit, Inhib, Astro, Oligo,
  OPC, Micro, EndoMural, Unattributed) keeps them from being adjacent.
* **Donor-level summaries** first average spot proportions within each
  donor x domain combination (one sample per donor), then summarize across
  donors. This weights donors equally rather than spots. Donor x domain
  combinations with fewer than 10 spots are left out of the plots (7 of 438,
  all in `SpD04-WM`, with 1-8 spots each). They remain in the saved csv
  together with `n_spots`.

## Merge QC

`01-merge_RESCUE_to_spe.R` writes every check to
`RESCUE_merge_qc_summary.csv`. *Hard* checks indicate misalignment and stop
the script if they fail; *diagnostic* checks are recorded only.
Per-sample spot coverage goes to `RESCUE_merge_qc_per_sample.csv`.

Results from the current run:

| Check | Type | Result |
| --- | --- | --- |
| Spots in `spe` / RESCUE | diagnostic | 279,806 / 279,806 |
| Duplicated RESCUE `spot_id` | hard | 0 |
| RESCUE spots matched to `spe` (by `key`) | diagnostic | 279,806 (100%) |
| `spe` spots with a RESCUE result | diagnostic | 279,806 (100%); every one of the 63 samples is 100% matched |
| Samples in only one source | diagnostic | none (63 in both) |
| Agreement of `sample_id`, `array_row`, `array_col`, `pxl_row_in_fullres`, `pxl_col_in_fullres` | hard | 100% identical |
| RESCUE `domain` vs `spe$PRECAST_07` | hard | 100% identical |
| No `NA` proportions; all within [0, 1] | hard | pass |
| Max deviation of per-spot proportion sum from 1 | hard | 2.8e-15 |
| RESCUE `nUMI` vs `spe$sum_umi` | diagnostic | Pearson r = 0.999999. Exact match in only 5.2% of spots: RESCUE `nUMI` is never larger (median ratio 0.9988, min 0.9797), consistent with RESCUE counting only the genes it used |
| `unattributed_ratio` vs `Unattributed Expression` | diagnostic | Not identical (r = 0.990; difference median -0.013, range -0.24 to 0.19). These appear to be two different RESCUE quantities, so both are kept |

**Conclusion:** the merge is a complete one-to-one match. Every spot in
both sources is paired with its counterpart, and spatial coordinates and
domain labels agree for every matched spot.

## Outputs

* `processed-data/rds/18_rescue/fnl_spe_kept_spots_only_w_RESCUE.rds`:
  finalized `spe` with RESCUE columns added to `colData`
* `processed-data/rds/18_rescue/RESCUE_prop_df.rds`: slim data frame
  (`key`, `sample_id`, `brnum`, `dx`, `spd_label`, `RESCUE_*`)
* `processed-data/rds/18_rescue/RESCUE_merge_qc_summary.csv`
* `processed-data/rds/18_rescue/RESCUE_merge_qc_per_sample.csv`
* `processed-data/rds/18_rescue/RESCUE_donor_spd_mean_prop.csv`: donor x
  domain x cell type mean proportions, with `n_spots`
* `plots/18_rescue/spot_plot_RESCUE_broad_class_rep_samples.pdf` (1 page
  per broad class)
* `plots/18_rescue/spot_plot_RESCUE_celltype_rep_samples.pdf` (1 page per
  cell type)
* `plots/18_rescue/RESCUE_broad_class_by_spd_rep_samples{,_excl_unattributed}.pdf`
* `plots/18_rescue/RESCUE_celltype_by_spd_heatmap_rep_samples.pdf`
* `plots/18_rescue/RESCUE_broad_class_by_spd{,_excl_unattributed}.pdf`:
  overall composition per spatial domain, diagnoses pooled
* `plots/18_rescue/RESCUE_broad_class_by_spd_dx{,_excl_unattributed}.pdf`
* `plots/18_rescue/RESCUE_broad_class_by_donor_spd{,_excl_unattributed}.pdf`
* `plots/18_rescue/RESCUE_celltype_by_spd_heatmap_dx.pdf`
* `plots/18_rescue/RESCUE_celltype_by_spd_heatmap_dx_row_scaled.pdf`
* `plots/18_rescue/RESCUE_donor_prop_by_spd_dx_boxplot.pdf` (page 1 =
  broad classes; then 1 page per broad class with >1 fine type)
