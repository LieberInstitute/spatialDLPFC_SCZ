# SCZ PTN/FGF mediation screening

This analysis adapts `LFF_spatial_ERC/code/22_Mediation` to matched donor–SpD
pseudobulks in spatialDLPFC_SCZ. It is an exploratory expression-association
screen, not a causal mediation estimator or functional validation experiment.

## Prespecified analyses

| Screen | Source mediator | Target | Observations | Donors |
|---|---|---|---:|---:|
| `ptn_vasc_neuropil` | vascular PTN | neuropil | 272 | 61 |
| `ptn_vasc_neun` | vascular PTN | neuronal | 220 | 58 |
| `fgf1_neuropil_vasc` | neuropil FGF1 | vascular | 272 | 61 |
| `fgf2_neuropil_vasc` | neuropil FGF2 | vascular | 272 | 61 |
| `fgf1_neun_vasc` | neuronal FGF1 | vascular | 220 | 58 |

All primary screens use SpD02/03/05/06/07. Neuronal FGF2 is excluded. Counts are
TMM-normalized, filtered with the matched baseline design, and fitted using
`edgeR::voomLmFit(adaptive.span=TRUE, sample.weights=TRUE)` with donor blocking.
The contrast is SCZ minus NTC; the primary covariates are age, sex, and SpD.
The source mediator is standardized normalized logCPM, one SD per matched set.

See [DESIGN.md](DESIGN.md) for decisions and caveats,
[IMPLEMENTATION.md](IMPLEMENTATION.md) for data and code contracts, and
[METHODS.md](METHODS.md) for draft manuscript methods.
[VALIDATION.md](VALIDATION.md) records the executed checks, primary findings, and
completion status.
[PLAN.md](PLAN.md) records scope and implementation phases;
[MODEL_SPECIFICATION.md](MODEL_SPECIFICATION.md) defines the formulas, contrasts,
variance models, and testing families;
[WORK_LOG.md](WORK_LOG.md) records implementation and execution;
[mediation-analysis-decisions-notes.md](mediation-analysis-decisions-notes.md)
records the analysis decisions and methodological references.
The [report archive](reports/2026-10-01/README.md) contains the completed run
tables, figures, logs, software provenance, and checksums within this folder.

The completed October 1, 2026 run tested 40,329 primary mediator–outcome pairs
and produced zero qualifying hits; all 17 required robustness analyses also
produced zero hits. All five matched source mediator diagnosis gates failed.
See [the run report](reports/2026-10-01/report/RESULTS.md)
and [the validation record](VALIDATION.md) for results and interpretation.

## Run on JHPCE

From the SCZ project root:

```bash
bash code/analysis/17_mediation/submit.sh
```

This submits a 4 CPU / 16 GB screening job. After it succeeds, a 4 CPU / 16 GB
sensitivity job and a 1 CPU / 64 GB raw-spot job run independently. A final
1 CPU / 8 GB report job waits for both. Analysis jobs initially request four
hours; reporting requests one hour. They use `conda_R/4.5`; requests are starting
resource budgets, not runtime guarantees. Existing inputs are read-only.
Only this analysis directory needs to be transferred to the remote repository.

Individual stages, executed within a compute allocation:

```bash
module load conda_R/4.5
Rscript --vanilla code/analysis/17_mediation/tests/test_core.R
bash code/analysis/17_mediation/tests/test_hit_paths.sh
Rscript --vanilla code/analysis/17_mediation/run_mediation.R --stage audit
Rscript --vanilla code/analysis/17_mediation/run_mediation.R --stage historical
Rscript --vanilla code/analysis/17_mediation/run_mediation.R --stage screen --workers 4
Rscript --vanilla code/analysis/17_mediation/run_mediation.R --stage sensitivity
Rscript --vanilla code/analysis/17_mediation/run_mediation.R --stage overlap
Rscript --vanilla code/analysis/17_mediation/run_mediation.R --stage report
python3 code/analysis/17_mediation/tests/verify_outputs.py --outdir processed-data/17_mediation --require-robustness
```

`--stage primary` runs audit, historical reconciliation, screening, sensitivities,
and reporting in one R session. `--stage all` also runs overlap and requires the
larger allocation. Default stage
is the read-only-data audit, which writes only analysis outputs. Other options:
`--project-root`, `--config`, `--outdir`, `--plots`, `--workers`, `--screens` and
`--audit-only`. `--screens` accepts comma-separated IDs. Use a separate output
folder for subset smoke runs; pooled five-screen FDR is unavailable for a subset.

R packages: SpatialExperiment, edgeR, limma, data.table, digest, and Matrix.
Tests are deterministic and require the same environment. The separate synthetic
hit-path fixture exercises whole-donor influence, between/within-donor models,
joint FGF models, raw count reconstruction, removal of shared spots from both
contexts, insufficient-spot handling, rejection of inconsistent raw counts, and
zero-hit control flow. Its forced triggers are confined to temporary
synthetic data and are never treated as scientific findings. The final report job independently verifies exported sample/gene identifiers,
shared baselines, per-screen/pooled BH corrections, and complete robustness
tables with current signatures and full gene universes using Python standard
libraries. No local package
installation is needed when working through JHPCE.

## Outputs

Outputs default to `processed-data/17_mediation`; figures to `plots/17_mediation`.

- `audit/`: input inventory, sample matching, exclusions, target gene universes.
- `historical/`: historical model reproduction and numerical reconciliation.
- `primary/all_pairs.tsv.gz`: complete results, including failed screening gates.
- `primary/screening_hits.tsv`: pairs satisfying the prespecified screening rule.
- `primary/mediator_gates.tsv`: source diagnosis effects, genome-wide FDR and
  supplementary adjustment across the five nominated candidate tests.
- `sensitivity/`: slide/RIN, matched logcounts, fixed weights, and neuronal
  SpD07 exclusion analyses, each fitted over its full retained gene universe.
- `diagnostics/`: hit-triggered leave-donor-out, between/within-donor, and
  neuropil FGF1/FGF2 joint-mediator results.
- `overlap/`: raw-to-pseudobulk count reconciliation, actual shared-spot counts,
  diagnosis summaries, and hit-triggered
  reaggregation with shared spots removed from both contexts.
- `report/RESULTS.md`, `screen_summary.tsv`, and hit sensitivity comparisons.
- `provenance/`: signatures, input/code fingerprints, configuration, versions.
- `cache/`: exact model-signature caches; `primary/*.rds`: reproducible fit objects.
- `logs/`: submitted job IDs and scheduler output.

A failed fit is recorded as failed, not as zero hits. A successful run with zero
hits is a valid scientific result. Check both sensitivity and overlap status
files before describing robustness as established. A report can be generated
before overlap finishes; it explicitly marks that stage as unverified.

## Screening rule

A pair must pass all five gates:

1. Historical target diagnosis nominal p < 0.05.
2. Matched voom baseline target diagnosis nominal p < 0.05.
3. Matched source mediator diagnosis nominal p < 0.05.
4. Conditional mediator–outcome BH FDR < 0.05 across all target genes in that screen.
5. Mediator-adjusted target diagnosis p >= 0.05.

Higher-priority hits also require same-direction coefficient shrinkage and
mediator–outcome BH FDR < 0.05 pooled across all five screens. This label does not
imply calibrated mediation-discovery FDR. Nominal mediator–outcome p-values are
provided but do not qualify a hit.

None of the four original mediator nominations passes genome-wide FDR 0.05.
Historical nominal p-values are PTN vascular .0131, FGF1 neuropil .0403,
FGF2 neuropil .0331, and FGF1 neuronal .0472. The screen estimates conditional associations; relative coefficient shrinkage
is not a causal proportion.
