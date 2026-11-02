# SCZ-DEG analysis in vasculature (Claudin-5+) spots

Differential expression between SCZ and NTC among `vasc_pos` spots, across **all
seven spatial domains** (no GM/WM subsetting), followed by GO BP GSEA.

`vasc_pos` spots have `spg_PClaudin5` between 0.05 and 0.20, as defined in
`code/analysis/01_build_spe/20_create_final_spe.R`.

## Scripts

Run from the repository root, in this order. CSV and `.rds` results are in
`processed-data/20_dx_deg_spg/dx_deg_spg_vasc/`, and PDFs are in
`plots/20_dx_deg_spg/dx_deg_spg_vasc/`. See `../README.md` for the overall layout.

| Order | Script | Purpose | Output |
|---|---|---|---|
| 1 | `00-create_pb_data.R` | Pseudobulk `vasc_pos` spots by `fnl_spd` and `sample_id` (`min_ncells = 10`) | `processed-data/rds/PB_dx_spg/pseudo_vasc_pos_donor_spd_fnl.rds` |
| 2 | `01-vasc_dx-DEG_allSpD.R` | DE test and exploratory checks | `vasc-dx_DEG-allSpD.csv`, `vasc-dx_DEG-allSpD_fdr_05.csv` |
| 3 | `GSEA_vasc_dx-DEG_allSpD.R` | GO BP GSEA, genes ranked by `t_stat_scz` | `gsea_vasc-dx_DEG-allSpD.csv`, `gsea_vasc-dx_DEG-allSpD.rds` |
| 4 | `treemap_GSEA_vasc_dx-DEG_allSpD.R` | Reduce GO terms with rrvgo and draw treemaps | `treemap_GO_vasc-dx_DEG-allSpD{,_up,_down}.pdf`, `rrvgo_reduced_terms_vasc-dx_DEG-allSpD.csv` |

## DE model

The model follows the `SPD` branch of
`code/analysis/99_archived/SPG_analysis/SPG_DE_analysis.R`:

- `registration_model()` with covariates `age` and `sex`
- `registration_block_cor()` with `sample_id` as the block
- `registration_stats_enrichment()` with covariates `age`, `sex` and `PRECAST_07`

## Results

- Input: 352 pseudobulks from 61 samples, 5,215 genes. 64 sample-domain
  pseudobulks were dropped for having fewer than 10 spots.
- DEGs: 11 at FDR < 0.05, 19 at FDR < 0.1.
- FDR < 0.05 genes: RPS26, TF, CEBPD, A2M, KIAA0930, DDIT4, MTRNR2L1, ALDH1A1,
  RASD1, GSTM5, ENPP2.
- Overlap with the 172 dxDEGs (`10_dx_deg_adjust_spd/172_prelim_fdr010.csv`):
  9 of the 11, all except CEBPD and KIAA0930.
- GSEA: 76 GO BP terms at FDR < 0.05, 60 negatively and 16 positively enriched.
  The top negative terms relate to mitochondrial electron transport and oxidative
  phosphorylation. The top positive terms are responses to copper and zinc ions.

These results reproduce `processed-data/spg_pb_de/test_SPD_pseudo_vasc_pos.csv`
exactly for all 5,215 genes (t-statistics and FDR).

## Notes

- `pseudo_vasc_pos_donor_spd.rds` in the same folder is the 2024 object, still
  read by `15_eqtl_coloc` and `code_for_isee_app`. It has the same genes, samples
  and counts as the `_fnl` object, but lacks some `colData` (e.g. `slide_id`,
  `fnl_spd`).
- The neun, neuropil and pnn folders test GM domains only, so their results are
  not directly comparable with this folder's.
- The final volcano plot in `01-vasc_dx-DEG_allSpD.R` is drawn but not saved.
