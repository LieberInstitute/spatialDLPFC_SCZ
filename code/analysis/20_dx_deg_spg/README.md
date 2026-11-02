# SCZ-DEG analysis within SPG-defined spots

Differential expression between SCZ and NTC, run separately on the spots positive
for each SPG channel, followed by GO BP GSEA.

| Folder | Spots | Spatial domains tested |
|---|---|---|
| `dx_deg_spg_neun` | `neun_pos` | GM only (spd02, 03, 05, 06, 07) |
| `dx_deg_spg_neuropil` | `neuropil_pos` | GM only (spd02, 03, 05, 06, 07) |
| `dx_deg_spg_pnn` | `pnn_pos` | GM only (spd02, 03, 05, 06, 07) |
| `dx_deg_spg_vasc` | `vasc_pos` | All seven domains |

The SPG calls are made in `code/analysis/01_build_spe/20_create_final_spe.R`.

## Where things live

Each channel has a folder of the same name in three locations. Scripts read and
write these paths with `here()`, so run them from the repository root.

| Content | Location |
|---|---|
| Code | `code/analysis/20_dx_deg_spg/<folder>/` |
| CSV results and GSEA `.rds` objects | `processed-data/20_dx_deg_spg/<folder>/` |
| PDF plots | `plots/20_dx_deg_spg/<folder>/` |
| Pseudobulk objects | `processed-data/rds/PB_dx_spg/` |

## Workflow

`<spg>` is the channel name and `<set>` is `GM`, or `allSpD` for vasculature.

| Order | Script | Purpose | Output |
|---|---|---|---|
| 1 | `00-create_pb_data.R` | Pseudobulk positive spots by `fnl_spd` and `sample_id` (`min_ncells = 10`) | `pseudo_<spg>_pos_donor_spd.rds` (`_fnl` suffix for vasculature) |
| 2 | `01-<spg>_dx-DEG_<set>.R` | DE test and exploratory checks | `<spg>-dx_DEG-<set>.csv`, `<spg>-dx_DEG-<set>_fdr_05.csv` |
| 3 | `GSEA_<spg>_dx-DEG_<set>.R` | GO BP GSEA, genes ranked by `t_stat_scz` | `gsea_<spg>-dx_DEG-<set>.csv`, `.rds` |
| 4 | `treemap_GSEA_<spg>_dx-DEG_<set>.R` | Reduce GO terms with rrvgo and draw treemaps | `treemap_GO_<spg>-dx_DEG-<set>{,_up,_down}.pdf`, `rrvgo_reduced_terms_<spg>-dx_DEG-<set>.csv` |

`dx_deg_spg_neuropil` has three more scripts:

| Script | Output |
|---|---|
| `02_plot_volcano-neuropil_dx-DEG_GM.R` | `volcano_plot_neuropil-dx_DEG-GM.pdf` |
| `volcano_plot_neuropil_dx_res.R` | `neuropil_dx_DEG_GM.pdf` |
| `Enrichement_GO_neuropil_dx.R` | `GO_neuropil_GW_up_gene.csv`, `GO_neuropil_GW_down_gene.csv`, `GSEA_neuropil_GW_down_gene.csv` |

## DE model

All four channels use the same model on the pseudobulked data:

- `registration_model()` with covariates `age` and `sex`
- `registration_block_cor()` with `sample_id` as the block
- `registration_stats_enrichment()` with covariates `age`, `sex` and `PRECAST_07`

## Results

| Channel | Genes tested | FDR < 0.05 | FDR < 0.1 | GSEA terms (FDR < 0.05) | NES > 0 | NES < 0 |
|---|---:|---:|---:|---:|---:|---:|
| NeuN | 16,768 | 16 | 90 | 124 | 15 | 109 |
| Neuropil | 11,834 | 145 | 352 | 400 | 243 | 157 |
| PNN | 10,879 | 0 | 0 | 55 | 4 | 51 |
| Vasculature | 5,215 | 11 | 19 | 76 | 16 | 60 |

Vasculature uses a different domain set from the other channels, so its counts
are not directly comparable with theirs.

## Other code that reads these results

- `code/analysis/15_eqtl_coloc/table_utils.R` (neun, neuropil and pnn)
- `code/analysis/prs_deg_spg_neuropil/03_cor_plot_neuropil_PRSvsDx_DEG.R`
- `code/image_processing/revisions/dx_deg_neuropil.R`

`table_utils.R` and the eQTL `.Rmd` files read vasculature results from
`processed-data/spg_pb_de/test_SPD_pseudo_vasc_pos.csv`. That file is numerically
identical to `vasc-dx_DEG-allSpD.csv` here.
