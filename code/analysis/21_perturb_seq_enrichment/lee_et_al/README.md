# Perturb-seq enrichment: Lee et al. (2026)

Tests whether the gene expression changes caused by knocking out SCZ risk genes in
mouse neocortex match our SCZ-DEG results from postmortem DLPFC.

> Status: implemented. See [Scripts](#scripts) and [Results](#results).

## Source study

Lee J, Min H, ... Won H. *Perturbation of genes linked to common schizophrenia
risk variants identifies cilia programs.* bioRxiv 2026.
[doi:10.64898/2026.06.09.731172](https://doi.org/10.64898/2026.06.09.731172)
([PMC13278005](https://pmc.ncbi.nlm.nih.gov/articles/PMC13278005/))

- **Design:** in vivo CROP-seq. Each of 12 SCZ risk genes was knocked out in the
  P14 mouse neocortex, read out by scRNA-seq.
- **Knocked-out genes:** Banf1, Dpf2, Eif1ad, Fxr1, Matr3, Mef2c, Nsd3, Paip2,
  Rbm26, Rnaseh2c, Scyl1, Srebf1.
- **Code:** <https://github.com/thewonlab/SCZ_CROPseq>. Only the cNMF analyses are
  shared. Their postmortem enrichment test is
  `cNMF/scripts/06_scz_postmortem_enrichment.R`.
- **Use of our data in their paper:** Fig 6C runs fgsea of their Matr3 and Scyl1
  DEGs against our PRS and per-domain SCZ t-statistics.

## Input data

`processed-data/perturb-seq-degs/2_won/R01_12_perts_full_glm_result.rda` holds the
object `full_glm_result`. `Won_README_SHK.docx` in the same folder describes it.

- 1.74M rows. Each row is one cell-type grouping, one perturbation and one tested
  gene, with both significant and non-significant genes.
- Negative binomial GLM:
  `target ~ gRNA + nCount_RNA + nFeature_RNA + nCount_crispr + nFeature_crispr + batch_num + percent.mt`
- 8 cell-type groupings: `Allcells`, `Neuron`, `Glutamatergic`, `GABAergic`,
  `Astro`, `Oligo`, `Microglia`, `Vascular`. `Allcells` and `Neuron` contain the
  other groupings.

### Known caveats

- **Not the paper's final DEG list.** `fdr < 0.05` gives 3,307 unique `Allcells`
  DEGs, while the paper reports 3,031. The paper's 1,000-permutation filter is
  not in the file.
- **Direction comes from `glm_coef`, not `raw_logfc`.** Perturbed cells have
  higher sequencing depth, so the raw fold change disagrees in sign with
  `glm_coef` for about 94% of the down DEGs.
- **Up DEGs are mostly near-zero genes.** Significant up genes have a median
  control mean of about 1e-4 UMI per cell. 326 rows have an absolute coefficient
  above 10, mostly in `Vascular` and `Oligo`.
- **Cell-type-level results are extra.** The paper only reports `Allcells`. The
  cell-type-level results are mostly upregulated genes.
- **Sanity check passes.** Each knocked-out gene is significantly down in its own
  perturbation in `Allcells`.

## Decisions

| Item | Decision |
|---|---|
| Their DEG definition | `fdr < 0.05`, split into up and down by the sign of `glm_coef` |
| Cell-type groupings | All 8. `Allcells` and `Neuron` are shown in a separate panel because they contain the others |
| Our gene sets | SCZ-DEG and PRS results only. Domain markers are deferred |
| Main test | fgsea |
| Orthologs | `babelgene` (the database behind orthogene's babelgene method; orthogene itself is not installed), `top = FALSE`, strict one-to-one only, mapped to human Ensembl ID. This gives 15,573 of 25,237 tested mouse genes |
| Our DEG cutoff (secondary analysis) | FDR < 0.10 |

## Main analysis: their gene sets against our rankings

**Gene sets.** 12 perturbations × 8 cell-type groupings × 2 directions = 192
sets, before the size filter.

**Our rankings**, each ranked by its moderated t:

| Family | File | Rank by | Number of rankings |
|---|---|---|---|
| Domain-adjusted SCZ | `processed-data/rds/10_dx_deg_adjust_spd/dx-deg_PRECAST07.csv` | `t_stat_scz` | 1 |
| Per-domain SCZ | `processed-data/rds/11_dx_deg_interaction/layer_restricted_degs_all_spds.csv`, split by `PRECAST_spd` | `t` | 7 |
| SPG-specific SCZ | `processed-data/20_dx_deg_spg/dx_deg_spg_{neun,neuropil,pnn}/<spg>-dx_DEG-GM.csv` and `dx_deg_spg_vasc/vasc-dx_DEG-allSpD.csv` | `t_stat_scz` | 4 |
| PRS | `processed-data/rds/14_prs_deg/norm_PRS_DEG_test_res_PRECAST07_donor_spd.csv` | `t` | 1 |

That makes 13 rankings.

**Settings**

1. **Gene universe.** For each perturbation and cell-type grouping, our ranking
   is limited to genes that were also tested in that mouse model, after ortholog
   mapping. This means up to 96 universes × 13 rankings.
2. **Near-zero up genes.** Matching to our tested genes should remove many of
   them, because our data already filtered out low-expression genes. Report how
   many genes each set keeps, and only add an expression filter if the up sets
   still look dominated by them.
3. **Set sizes.** `minSize = 15`, and `maxSize` about 1,000, since the largest
   `Allcells` down sets have around 780 genes.
4. **Multiple testing.** fgsea's `padj` corrects within one ranking. Also apply BH
   across each family's full grid of rankings × gene sets.
5. **Reading the results.** A perturbation down set with a negative NES means
   the knockout pushes expression in the same direction as SCZ or higher PRS.

## Check against their figure

Before the full grid, run `Allcells` Matr3 and Scyl1 against the per-domain and
PRS rankings. Expected: enrichment among genes that are down across L2–6, plus
Scyl1 down in L1 and up in WM, as in their Fig 6C. Our gene sets are a slight
superset of theirs, so close agreement is enough.

## Secondary: our DEG sets against their rankings

Rank their genes per perturbation and cell-type grouping by
`sign(glm_coef) × -log10(glm_p)`, then test our up and down DEG sets against it.
Many of our sets are small and will fall under `minSize`, so this is only a
check.

## Output

- One NES heatmap per family of our rankings: their gene sets as rows, grouped
  by cell type, and our rankings as columns, with asterisks for FDR.
- Leading-edge genes for the significant results.
- A table of how many genes each set keeps after ortholog mapping and universe
  matching.

Cell types that line up with our data, for interpretation only:

- `Neuron` and `Glutamatergic` with NeuN
- `GABAergic` with PNN
- `Vascular` with vasculature
- `Oligo` with WM

## Layout

Code lives in this folder, `code/analysis/21_perturb_seq_enrichment/lee_et_al/`.
Outputs follow the repo pattern: `processed-data/21_perturb_seq_enrichment/lee_et_al/`
and `plots/21_perturb_seq_enrichment/lee_et_al/`. The other two perturbation
papers will get sibling folders.

## Scripts

Run from the repository root, in order, or with `bash run_all.sh` from this
folder. Each step runs locally in under 5 minutes. Logs go to `logs/`.

| Step | Script | Output (`processed-data/...` or `plots/...`) |
|---|---|---|
| 00 | `00-prep_lee_gene_sets.R`: ortholog mapping, gene sets, universes | `lee_ortholog_1to1.csv`, `lee_ortholog_status.csv`, `lee_deg_mapping.csv`, `lee_glm_human.rds`, `lee_gene_sets.rds`, `lee_gene_set_size.csv` |
| 01 | `01-prep_our_rankings.R`: our 13 t-statistic rankings | `our_rankings.rds` |
| 02 | `02-fgsea_lee_sets_vs_our_rankings.R`: main analysis | `fgsea_lee_sets_vs_our_rankings.{csv,rds}` |
| 03 | `03-plot_fgsea_lee_sets_vs_our_rankings.R` | `fig6C_check_Allcells_Matr3_Scyl1.pdf`, `heatmap_NES_{domain_adjusted_PRS,per_domain,spg}.pdf` |
| 04 | `04-fgsea_our_sets_vs_lee_rankings.R`: secondary analysis | `fgsea_our_sets_vs_lee_rankings.{csv,rds}`, `our_deg_set_size.csv`, `our_deg_sets.csv` |
| 05 | `05-plot_fgsea_our_sets_vs_lee_rankings.R` | `heatmap_NES_our_sets_vs_lee_rankings.pdf` |
| 06 | `06-plot_ortholog_matching.R`: how genes are kept or dropped by ortholog and universe matching | `ortholog_matching_summary.pdf` |
| 07 | `07-fig6C_check_all_DEGs.R`: Fig 6C check with up and down DEGs pooled, tested signed (t) and undirected (\|t\|, `scoreType = "pos"`) | `fig6C_check_all_DEGs.csv`, `fig6C_check_Allcells_Matr3_Scyl1_all_DEGs.pdf` |
| 08 | `08-fgsea_lee_pooled_sets_vs_our_rankings.R`: main analysis with up and down DEGs pooled, signed and undirected | `fgsea_lee_pooled_sets_vs_our_rankings.{csv,rds}` |
| 09 | `09-plot_fgsea_lee_pooled_sets_vs_our_rankings.R` | `heatmap_NES_pooled_{domain_adjusted_PRS,per_domain,spg}.pdf` |
| 10 | `10-fgsea_our_pooled_sets_vs_lee_rankings.R`: secondary analysis with our up and down DEGs pooled, signed and undirected | `fgsea_our_pooled_sets_vs_lee_rankings.{csv,rds}`, `our_pooled_deg_set_size.csv` |
| 11 | `11-plot_fgsea_our_pooled_sets_vs_lee_rankings.R` | `heatmap_NES_our_pooled_sets_vs_lee_rankings.pdf` |
| 10, 11 (`p05`) | Same as 10 and 11 with our DEGs at nominal p < 0.05 and maxSize 3000: `Rscript <script> p05` | `fgsea_our_pooled_sets_vs_lee_rankings_nomp05.{csv,rds}`, `our_pooled_deg_sets_nomp05.csv`, `heatmap_NES_our_pooled_sets_vs_lee_rankings_nomp05.pdf` |

**Heatmap layout.** In the main-direction heatmaps (03, 07, 09) rows are their
gene sets (cell-type grouping × perturbation; down and up side by side in
03) and columns are our rankings. In the reverse heatmaps (05, 11) rows are
our DEG sets and columns are their cell-type grouping × perturbation. Every
heatmap has stacked bars on the rows and columns, split into genes with and
without a one-to-one ortholog (`utils_heatmap_bars.R`):

- Their DEG sets (rows in 03/07/09): DEGs with a one-to-one human ortholog vs not.
- Our rankings (columns in 03/07/09): ranked genes with a one-to-one match to a
  mouse gene tested in CROP-seq vs not.
- Our DEG sets (rows in 05/11): same definition as our rankings.
- Their rankings (columns in 05/11): genes tested in that model, with a
  one-to-one human ortholog vs not.

The bars show set sizes before universe matching. Figures label our rankings "Layer-adjusted" and "Layer-restricted" (data files keep
"Domain-adjusted" / per_domain). The pooled figures (07, 09,
11) show the signed test only; the undirected results are kept in the csv/rds
files.

`utils_heatmap_bars.R` is sourced by the plotting scripts and is not run on
its own.

Pooled steps (07–11) add a `version` column. `signed` ranks by t (or their
signed score) with `scoreType = "std"`. `undirected` ranks by |t| (or
-log10 p) with `scoreType = "pos"`, so its NES is always positive. BH is
applied within each version.

Result columns: `padj_ranking` is BH within one of our rankings; `padj_family`
is BH across a whole family and is the one used in the plots. `tested` means
the set passed the size filter. `retried` means fgsea was rerun with
`nPermSimple = 10000` after an NA p-value.

## Results

First pass, 2026-11-12.

**Gene-set sizes.** Ortholog mapping keeps about 85% of `Allcells` down DEGs
but only 25–40% of up DEGs. Most unmapped up genes are sparse `Gm` and `Rik`
genes. Of the 25,237 tested mouse genes, 15,573 (62%) have a one-to-one
human ortholog. The other genes have no human ortholog (8,500), map to
several human genes (690), or share a human gene with another mouse gene
(474). 6,909 of the 8,500 genes without an ortholog are `Gm` or `Rik` genes.
Between 4,820 (Vasc) and 12,782 (NeuN) of the matched genes are also tested in
our rankings. `ortholog_matching_summary.pdf` shows each step. After universe matching and `minSize = 15`, 1,126 of 2,444 tests remain.
Most cell-type-level down sets in `Microglia`, `Oligo` and `Vascular` are too
small to test.

**Check against Fig 6C: reproduced.** `Allcells` Matr3 and Scyl1 down sets have
NES −2.1 to −3.1 across L2–6 and about −2.6 for PRS. Scyl1 down is negative in
L1 and positive in WM (NES 1.4). On the first pass fgsea failed to give
p-values for Matr3, Scyl1 and Fxr1 in `SpD01-WMtz`, matching the grey boxes in
their figure. With the retry, Scyl1 and Fxr1 are significant there and Matr3
is not.

**Main analysis.** Most significant results have negative NES, meaning the
perturbation down sets sit among genes that are down in SCZ or with higher
PRS:

| Family | Tests | padj_family < 0.05 | Negative NES | Positive NES |
|---|---:|---:|---:|---:|
| Domain-adjusted | 91 | 39 | 38 | 1 |
| Per-domain | 637 | 205 | 189 | 16 |
| SPG | 307 | 102 | 97 | 5 |
| PRS | 91 | 30 | 27 | 3 |

- All 12 `Allcells` down sets are significant with negative NES for
  domain-adjusted SCZ, PRS, NeuN, Neuropil and Vasc; PNN has 10 of 12. In the
  per-domain results all 12 are significant in L2/3, L3/4 and L5 and 11 in L6,
  with fewer in L1 (9) and WMtz (4). WM reverses (median NES +0.8, with 2
  significant negative).
- At the cell-type level, the clearest signal is in `Glutamatergic` and `Neuron`
  down sets. Glutamatergic Banf1 and Scyl1 are significant with negative NES
  in 12 of 13 rankings, and Mef2c and Nsd3 in 11.
- Up sets are mostly not significant (39 of 575 tests).

**Secondary analysis.** 10 of our 26 up/down DEG sets have at least 15 genes in
the universe: domain-adjusted up and down, PRS up and down, NeuN down,
Neuropil up and down, SpD04-WM up and down, and SpD01-WMtz up. Neuropil down
(266 genes) and PRS down (240 genes) are enriched among genes knocked down by
most perturbations in `Allcells` and `Neuron`. This matches the main analysis
from the other direction. One test (`SpD04-WM up` against `Microglia` Nsd3)
still has an NA p-value after the retry.

**Pooled up + down DEG sets (steps 07–11), 2026-11-18.** Pooling makes 96
sets instead of 192, and 1,002 tests pass the size filter in each version
(1,126 for the split sets). Results at `padj_family < 0.05`:

| Family | Tests | Signed: significant (negative / positive) | Undirected: significant |
|---|---:|---:|---:|
| Domain-adjusted | 80 | 36 (36 / 0) | 24 |
| Per-domain | 560 | 178 (168 / 10) | 117 |
| SPG | 282 | 90 (89 / 1) | 64 |
| PRS | 80 | 24 (24 / 0) | 27 |

- The signed pooled results closely follow the split down-set results,
  because down genes dominate most pooled sets. `Allcells` is significant
  for 129 of 156 signed tests and 101 of 156 undirected tests. The per-domain
  pattern is the same: strongest in L2/3–L6, weaker in L1 and WMtz, with
  WM median NES +0.9 under the signed test.
- The undirected test is significant less often. Only for PRS is it as
  sensitive as the signed test. The clearest cell-type-level sets are still
  `Glutamatergic` Nsd3, Banf1, Mef2c and Scyl1.
- Pooling lets more cell-type-level sets reach `minSize`, so `Astro`,
  `Oligo` and `Microglia` now show a few signed hits (e.g., `Oligo` Eif1ad
  in 10 of 13 rankings).
- One signed test has no p-value even after the retry: `Allcells` Scyl1
  against SpD01-WMtz.
- **Reverse direction.** Pooled, 7 of our 13 rankings have enough DEGs to
  test (PNN, L3/4, L6, L5, L2/3 have 0–2 DEGs; L1 has 17, but fewer than 15
  are in any universe). Neuropil (352 genes) and PRS (354) are the strongest:
  67 and 57 of 96 signed tests significant, mostly negative, i.e., our DEGs
  concentrate among genes the knockouts push down.

**Reverse direction, our DEGs at nominal p < 0.05 (2026-11-18).** With the
looser cutoff every one of our 13 rankings has a testable pooled set (283–2,521
genes; maxSize raised to 3000). Signed test, `padj_family < 0.05`, out of 96
tests each: layer-adjusted 59, L1 24, L2/3 64, L3/4 58, L5 61, L6 32, WMtz 5,
WM 24, NeuN 61, Neuropil 69, PNN 29, Vasc 31, PRS 55; almost all negative.
All 12 `Allcells` perturbations are significant for layer-adjusted, L2/3,
L3/4, L5, NeuN, Neuropil and PRS. 13 tests have NA p-values after the retry.
