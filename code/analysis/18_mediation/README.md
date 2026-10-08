# PTN/FGF cross-microenvironment mediation screen

Baron & Kenny-style screen of whether SCZ-associated expression of a ligand in one
SPG microenvironment accounts for SCZ-associated expression in another, using the
same donor x SpD pseudobulks and DE model as the manuscript (Fig 6F LR pairs).
Adapted from `LFF_spatial_ERC/code/22_Mediation`. Findings: [REPORT.md](REPORT.md).

## Screens (`screens.tsv`)

| screen | mediator (source context) | outcomes (target context) |
|---|---|---|
| `ptn_vasc_neuropil` | PTN (vascular) | neuropil genes |
| `ptn_vasc_neun` | PTN (vascular) | neuronal genes |
| `fgf1_neuropil_vasc` | FGF1 (neuropil) | vascular genes |
| `fgf2_neuropil_vasc` | FGF2 (neuropil) | vascular genes |
| `fgf1_neun_vasc` | FGF1 (neuronal) | vascular genes |

## Model (`engine = "limma_logcounts"`, the manuscript SPG DE model)

- Observations: donor x SpD pseudobulks present in both contexts (>= 10 spots each),
  gray-matter SpD02/03/05/06/07. Stored `logcounts` (TMM log2-CPM, prior.count 1).
- a (Dx -> M): source logcounts `~ 0 + Dx + age + sex + SpD`, SCZ - NTC.
- c (Dx -> Y): target logcounts, same design. c' and b: same design `+ M`,
  where M is the source mediator logCPM standardized over matched observations.
- `limma::lmFit` with donor blocking; consensus correlation from
  `duplicateCorrelation` on `~ Dx + age + sex` (as `spatialLIBD::registration_block_cor`);
  one target correlation shared by the c and c'/b fits; default `eBayes`.
- B&K steps, as in the ERC framework:
  1. X -> Y: target gene is a manuscript microenvironment DEG (p < .05) and has
     matched-sample c p < .05 (c and c' come from the same samples).
  2. X -> M: accepted from the manuscript nomination of the mediator as a nominal
     DEG in its source context (PTN vascular p = .013, FGF1 neuropil .040,
     FGF2 neuropil .033, FGF1 neuronal .047). The matched-sample a-path is
     reported but does not gate hits (`require_mediator_gate = FALSE`).
  3. Y ~ X + M: b BH FDR < .10 over all target genes in the screen, and c' p >= .05.
  `higher_priority` also requires |c'| < |c| (same sign) and b FDR < .10 pooled over screens.
- Historical stage re-fits the three manuscript microenvironment DE tables exactly.
- Overlap stage: SPG labels are not exclusive, so spots labeled in both contexts are
  removed from both pseudobulks, which are rebuilt from the raw spots and refit
  (`overlap_refit_all = TRUE`). The raw spots must first reproduce the stored
  pseudobulk counts exactly.
- `ptn_by_spd.R`: vascular PTN SCZ - NTC per SpD (all 7 SpDs, interaction model).
- `engine = "voom"` (`--engine voom`) and `--stage sensitivity` are retained in code
  but are not part of the reported analysis.

## Run

Inputs (read-only, paths in `config.R`): the three `PB_dx_spg/pseudo_*_pos_donor_spd.rds`
objects, the manuscript DE tables, `processed-data/ref/donor_meta.tsv.gz`, and for the
overlap stage `01_build_spe/fnl_spe_kept_spots_only.rds` (2.6 GB, ~25 GB RAM).
R >= 4.5 with SpatialExperiment, edgeR, limma, data.table, digest, Matrix.

```bash
# local (from the project root), ~25 min with 5 workers
SLURM_CPUS_PER_TASK=5 bash code/analysis/18_mediation/job.sh \
  "$PWD/code/analysis/18_mediation" "$PWD" "$PWD/processed-data/18_mediation/<run>"
# JHPCE: one 5 CPU / 64 GB job
bash code/analysis/18_mediation/submit.sh
```

`job.sh` runs the unit tests, `run_mediation.R --stage all`
(audit, historical, screen, overlap, report), `ptn_by_spd.R`, and
`tests/verify_outputs.py` (recomputes all BH corrections from exported p-values).
`tests/test_hit_paths.sh` exercises the hit-triggered code paths on synthetic data.

## Outputs (`--outdir`)

- `primary/all_pairs.tsv.gz`: every screen x target gene with a/c/c'/b estimates,
  gates, and `failed_gates`; `mediator_gates.tsv`: one a-path row per screen.
- `overlap/`: raw-to-pseudobulk reconciliation, shared-spot counts,
  `<screen>_results.tsv.gz` refits and `<screen>_hit_comparison.tsv`
  (primary vs shared-spot-removed for primary hits).
- `historical/`, `audit/`: manuscript-model reproduction, matched samples, exclusions.
- `report/RESULTS.md`, `report/screen_summary.tsv`, `ptn_vasc_by_spd.tsv`.
- `provenance/<signature>/`: input/code fingerprints, config, sessionInfo.
