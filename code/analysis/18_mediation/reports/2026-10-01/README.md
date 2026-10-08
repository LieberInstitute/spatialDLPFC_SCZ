# Completed run archive: 2026-10-01

This directory preserves the available completed-run outputs within the analysis
folder. Runtime output directories outside this folder were read without changes.
The archive is a review snapshot, not an input directory for new fits.

## Contents

- `report/RESULTS.md`, `report/screen_summary.tsv`: primary overview.
- `report/sensitivity_summary.tsv`: derived counts, source diagnosis tests, and
  hit totals for all 17 required sensitivities.
- `report/design_columns.tsv`: model-matrix columns and ranks reconstructed
  during packaging from archived primary samples using the unchanged design
  function; includes baseline/adjusted and slide/RIN models.
- `primary/`: all 40,329 result rows, individual screen tables, sample metadata,
  mediator gate table, zero-row hit table with schema, and fit status/cache keys.
- `sensitivity/`: all 17 full-gene result tables and completion status.
- `audit/`: available and matched samples, exclusions, model ranks, and gene sets.
- `historical/`: complete reference-model reproduction and reconciliation.
- `overlap/`: raw/pseudobulk reconciliation, shared-spot counts, summaries,
  input fingerprint, and no-hit reaggregation statuses.
- `figures/primary_diagnostics.pdf`: five primary diagnostic pages.
- `logs/`: retained execution/integration logs and review-time verification logs.
- `provenance/`: code/input fingerprints, configuration, screen manifest, package
  versions, stage invocations, and remote checkout revision.
- `archive_manifest.tsv`: 77 imported artifacts with original and archived
  SHA-256 values, sizes, project-relative source paths, and normalization flags.

The derived sensitivity summary, reconstructed design columns, this README,
and review-time verification logs
were added during packaging and are not imported artifacts in that manifest.

## Path normalization and unchanged results

Absolute project paths in four metadata files are replaced with `<PROJECT_ROOT>`:
`primary/status.tsv`, `overlap/raw_input_fingerprint.tsv`, and the code/input
fingerprint tables in `provenance/`. Account home paths, if present, are normalized
to `<HOME>`. This removes machine/account-specific attribution from the snapshot.
The manifest records the original hashes as well as the changed archive hashes.
Numeric test results, anonymized scientific sample identifiers, source data MD5
values, cache keys, and the executed model signature are unchanged.

The recorded model signature describes the original runtime path-dependent
provenance, not a newly computed signature over these normalized metadata files.
The archived fingerprints describe the code and paths used for the original run.
The analysis and runtime output folders were subsequently renamed from
`17_mediation` to `18_mediation`, and the entry point defaults were updated.
Archived metadata and hashes retain the original paths and fingerprints. Python
packaging tools and Markdown documentation are outside the R model fingerprint.

## Reproduction and validation

From the repository root, the exported statistical results can be checked without
raw data or R packages:

```bash
python3 code/analysis/18_mediation/tests/verify_outputs.py \
  --outdir code/analysis/18_mediation/reports/2026-10-01 --require-robustness
```

New full-data fits require the input files identified in `config.R` and the
recorded Bioconductor environment. Default runtime outputs remain in
`processed-data/18_mediation` and `plots/18_mediation`; this packaging operation
does not change pipeline output defaults. The [main README](../../README.md)
documents stage commands and scheduler use.

`export_report_snapshot.py --label LABEL` copies available outputs into a new
snapshot directory and writes its manifest. Existing snapshots are not
overwritten. Raw expression objects, large fitted-model RDS files, and caches
are not duplicated into this snapshot. The small provenance configuration RDS
is included. Runtime fit-object paths in `primary/status.tsv` point to objects
retained on JHPCE; they do not imply those objects are in this archive.
