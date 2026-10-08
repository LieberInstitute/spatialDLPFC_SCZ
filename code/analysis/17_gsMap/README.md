# 17 - gsMap: mapping SCZ genetic risk onto Visium spots

This directory runs [gsMap](https://github.com/JianYang-Lab/gsMap)
(Song et al., *Nature* 2025, "Spatially resolved mapping of cells associated
with human complex traits") on our 63 Visium DLPFC samples. gsMap tests each
spot for enrichment of GWAS heritability (stratified LDSC) in the genes that
are specifically expressed in that spot. Gene specificity is computed from the
spot's neighbourhood in a graph-neural-network latent space plus physical
space. Spot-level p-values are then aggregated per spatial domain by Cauchy
combination.

Questions addressed:

1. Which spatial domains (PRECAST-07 `spd_label`) are enriched for SCZ
   GWAS heritability (PGC3)?
2. Does the spatial distribution of SCZ-risk-associated spots differ
   between SCZ and NTC donors?

## Data Inputs

* `processed-data/rds/01_build_spe/fnl_spe_kept_spots_only.rds`: finalized
  Visium `spe` (63 donors, QC-passed spots, `spd_label` annotated), produced
  in [01_build_spe](/code/analysis/01_build_spe).
* `processed-data/ref/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz`:
  PGC3 SCZ GWAS, European ancestry, hg19 (Trubetskoy et al. 2022); the same
  file used in [15_eqtl_coloc](/code/analysis/15_eqtl_coloc).
* gsMap reference resources (`gsMap_resource.tar.gz` from the Yang lab):
  1000G EUR Phase 3 plink LD panel, HapMap3 SNP list, LDSC regression weights,
  and gencode v46lift37 GTF. All are hg19 and so match PGC3. Downloaded by
  `00-setup_gsMap_env.sh`.

## Design decisions

| Decision | Choice | Rationale |
| --- | --- | --- |
| Unit of analysis ("slice") | One h5ad per Visium capture area (`sample_id`), 63 total | gsMap builds spatial neighbourhoods within a slice; per-sample jobs also parallelize as SLURM arrays. |
| Cross-sample comparability | `gsmap create_slice_mean` over all 63 samples, passed as `--gM_slices` | gsMap's recommended approach for multiple slices/replicates: gene specificity scores share a common reference across samples, so donors can be compared. |
| Annotation | `spd_label` (PRECAST-07 spatial domains) | Supervises the GNN latent representation and defines the groups for the Cauchy combination. |
| Expression input | Raw UMI counts in `layers['count']` (`--data_layer count`) | gsMap normalizes and selects HVGs internally from counts. |
| Gene identifiers | Gene symbols (`rowData(spe)$gene_name`), first occurrence kept for duplicated symbols | gsMap links genes to SNPs by matching symbols to its GTF. |
| SNP-to-gene mapping | TSS ± 50 kb ("TSS only" strategy) | gsMap default. The enhancer-linked strategies (ABC/Roadmap) are possible later sensitivity analyses. |
| GWAS sample size | Per-SNP N_eff = 4 / (1/NCAS + 1/NCON), computed from case/control counts | The standard effective N for case-control LDSC, computed explicitly rather than relying on the PGC3 `NEFF` column. |
| Traits | `SCZ_PGC3` (primary). More traits are added as rows in `traits.tsv`. | Commented-out templates for BIP, MDD and a height negative control are included. Each needs its raw sumstats formatted in `02-format_sumstats.sh`. |

## Workflow

All shared paths and parameters are in `gsMap_config.sh`, which every step
sources. Steps marked "array" run one SLURM task per sample (`--array=1-63`),
where task *N* processes line *N* of `sample_list.txt`.

| Step | Script | Type | Description |
| --- | --- | --- | --- |
| 00 | `00-setup_gsMap_env.sh` | interactive, once | Create conda env (`gsMap==1.73.8`, Python 3.11), record `pip freeze`, download and check gsMap resources. |
| 01 | `01-export_spe_to_h5ad.R` / `.sh` | single job | Export `spe` to one h5ad per sample: `layers['count']`, `obsm['spatial']`, `obs` (`key`, `brnum`, `dx`, `age`, `sex`, `spd_label`), with `obs_names = key` and gene symbols as `var_names`. Writes `sample_list*.txt` (all / NTC / SCZ) and `h5ad_list.txt`. |
| 02 | `02-prep_PGC3_sumstats.py`, `02-format_sumstats.sh` | single job | Clean PGC3 (drop `##` header, compute N_eff), then `gsmap format_sumstats` (INFO ≥ 0.9, MAF ≥ 0.01) to produce `SCZ_PGC3.sumstats.gz`. |
| 03 | `03-create_slice_mean.sh` | single job | `gsmap create_slice_mean` over all 63 samples, producing `spe_slice_mean.parquet`. |
| 04 | `04-find_latent_representations.sh` | array | `gsmap run_find_latent_representations`: GNN-VAE latent embedding (`latent_GVAE`). |
| 05 | `05-latent_to_gene.sh` | array | `gsmap run_latent_to_gene` with `--gM_slices`: per-spot gene specificity scores (51 latent / 201 spatial neighbours). |
| 06 | `06-generate_ldscore.sh` | array | `gsmap run_generate_ldscore` for chr 1-22: spot-level stratified LD scores. **Most compute-heavy step.** |
| 07 | `07-spatial_ldsc.sh` | array | `gsmap run_spatial_ldsc` for each trait in `traits.tsv`: spot-level enrichment p-values. |
| 08 | `08-cauchy_combination_per_sample.sh` | array | `gsmap run_cauchy_combination`: one p-value per spatial domain per sample and trait. |
| 09 | `09-cauchy_combination_across_samples.sh` | single job | Cross-sample Cauchy combination per domain, for all samples, NTC only, and SCZ only. |
| 10 | `10-report_rep_samples.sh` | single job, optional | gsMap HTML reports for the representative samples Br8667 (`V13M06-342_D1`, NTC) and Br5973 (`V13M06-343_D1`, SCZ). |
| 11 | `11-collect_gsMap_results.R` / `.sh` | single job | Collect gsMap outputs into tidy `key`-indexed spot results plus sample × domain and cross-sample summaries. |
| 12 | `12-viz_gsMap_results.R` / `.sh` | single job | Spot plots (rep samples), cross-sample domain bar plot, sample × domain heatmap, and SCZ vs NTC donor-level test per domain. |

Dependencies: 01 → {03, 04}; 04 + 03 → 05 → 06; 06 + 02 → 07 → {08, 09};
08 → 10; 08 + 09 → 11 → 12.

### How to run (JHPCE)

```bash
cd code/analysis/17_gsMap

# once, on an interactive compute node
bash 00-setup_gsMap_env.sh

# recommended: first test one sample end-to-end to tune resources, e.g.
#   sbatch 01-export_spe_to_h5ad.sh; sbatch 02-format_sumstats.sh; sbatch 03-create_slice_mean.sh
#   then 04-08 with --array=1
#   (sbatch --array=1 04-find_latent_representations.sh, etc.)

# full pipeline: submits 01-12 with SLURM dependencies
bash run_all.sh

# resume from a given step (earlier outputs already exist)
bash run_all.sh 07
```

`run_all.sh` chains the per-sample array steps (04 → 05 → 06 → 07 → 08) with
`aftercorr`, so each sample advances to the next step as soon as its own
previous task succeeds. Logs are written to `logs/`.

### Adding a trait

1. Format its raw sumstats with `gsmap format_sumstats` into
   `processed-data/17_gsMap/sumstats/<TRAIT>.sumstats.gz` (add the call to
   `02-format_sumstats.sh`).
2. Add a `<TRAIT>\t<TRAIT>.sumstats.gz` row to `traits.tsv`.
3. `bash run_all.sh 07`: LD scores (06) are trait-independent and are reused.

## Outputs

Large intermediate files (git-ignored) under `processed-data/17_gsMap/`:

* `h5ad/{sample_id}.h5ad`, `h5ad/sample_list*.txt`, `h5ad/sample_dx.csv`
* `sumstats/SCZ_PGC3.sumstats.gz`
* `slice_mean/spe_slice_mean.parquet`
* `workdir/{sample_id}/{find_latent_representations, latent_to_gene, generate_ldscore, spatial_ldsc, cauchy_combination, report}/`
  * spot-level: `spatial_ldsc/{sample_id}_{trait}.csv.gz` (`spot` = `spe$key`, `p`, ...)
  * domain-level: `cauchy_combination/{sample_id}_{trait}.Cauchy.csv.gz` (`annotation`, `p_cauchy`, `p_median`)
* `cauchy_across_samples/{trait}_{all,ntc,scz}.Cauchy.csv.gz`

Summarized results under `processed-data/rds/17_gsMap/`:

* `gsMap_spot_res_df.rds`: per spot × trait `p`, `neg_log10_p` (join to `colData(spe)` by `key`)
* `gsMap_sample_spd_df.rds` / `gsMap_sample_spd_summary.csv`: per sample × domain × trait: `n_spots`, `mean_neg_log10_p`, `prop_p_lt_0.05`, `p_cauchy`
* `gsMap_across_sample_cauchy_df.rds` / `.csv`
* `gsMap_dx_test_by_spatial_domain.csv`: SCZ vs NTC per domain (`mean_neg_log10_p ~ dx + age + sex`, BH-FDR within trait)

Plots under `plots/17_gsMap/`:

* `spot_plot_{trait}_rep_samples.pdf`
* `across_sample_cauchy_by_spatial_domain.pdf`
* `sample_by_domain_cauchy_heatmap_{trait}.pdf`
* `donor_mean_neg_log10_p_by_dx_{trait}.pdf`

## Notes and caveats

* **Interpretation.** gsMap measures whether a spot's *specifically expressed
  genes* are enriched for GWAS heritability. It does not measure a donor's
  genetic risk. Differences between SCZ and NTC donors (step 12, panel 4)
  reflect differences in expression programs, not genotype, and should be
  treated as exploratory. For donor genotype-based analyses see
  [14_prs_deg](/code/analysis/14_prs_deg).
* **Cauchy p-values depend on spot counts.** A per-sample Cauchy p-value
  becomes more significant as a domain contains more spots. The donor-level
  test therefore uses the mean spot-level −log10(p) instead.
* **Gene symbol matching.** Our reference is GRCh38 (gencode v32,
  `refdata-gex-GRCh38-2020-A`), while gsMap's GTF is gencode v46lift37. Genes
  whose symbols were renamed between releases will not map to SNPs. Check how
  many genes match in the step 05/06 logs.
* **Resources.** The SBATCH `--mem`/`--time` values are initial estimates
  for ~5k spots per sample. Tune them after a single-sample test run; step 06
  (LD scores) is the bottleneck.
* **Sensitivity analyses (not yet implemented).**
  * Enhancer-linked SNP-to-gene mapping (`--enhancer_annotation_file`).
  * A negative-control trait (e.g. height).
  * Related psychiatric traits (BIP, MDD).
