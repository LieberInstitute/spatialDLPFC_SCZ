# Implementation and adaptation guide

## Files

- `config.R`: input paths, domains, thresholds, minimum sample requirements,
  and audited expected counts.
- `screens.tsv`: fixed scientific hypotheses, gene IDs, expected directions,
  and matched sample counts.
- `run_mediation.R`: CLI, provenance, package checks, and stage dispatch.
- `R/core.R`: loading, matching, design matrices, count preparation, fit caching,
  contrasts, result assembly, and screening classification.
- `R/stages.R`: audit, historical reconciliation, primary fits, sensitivities.
- `R/overlap.R`: shared-spot counts and hit-triggered nonoverlapping pseudobulks.
- `R/report.R`: full-family summaries, robustness comparisons, and PDFs.
- `tests/test_core.R`: alignment, gates, multiplicity, design, and synthetic
  signal/negative-control checks.
- `tests/test_hit_paths.sh`: temporary synthetic integration fixtures for
  hit-triggered diagnostics and shared-spot rebuilding, including zero-hit behavior.
- `submit.sh`, `job.sh`: screening, parallel sensitivity/spot-audit jobs, and a
  final report with failure dependencies. Pass the analysis directory explicitly
  when submitting the standalone test script with sbatch, because Slurm changes
  the executing script location.

## Input contracts

A context object supplies `counts`, `logcounts`, `genes`, and `samples`.
Counts/logcounts are genes × matched observation matrices. Stable, unique,
version-stripped Ensembl IDs are used; ambiguous duplicated IDs are rejected.
Genes have `gene_id` and `gene_name`. Samples have `donor`, `section`, `Dx`,
`age`, `sex`, `SpD`, `ncells`, `slide_id`, `rin`, and `key = donor|SpD`.
Diagnosis labels normalize to NTC and SCZ. No expression or metadata imputation
is performed; vascular slide/RIN columns are recovered from the donor metadata
with cross-checks against paired context metadata.

Historical results require `ensembl`, `p_value_scz`, `fdr_scz`, and `logFC_scz`.
Historical nominal selection remains distinct from the matched voom baseline.
A screen specifies source/target context and one mediator gene. Candidates are
never inferred from a new p-value scan.

For a new dataset, explicitly revise the schemas, audited count expectations,
screen manifest, donor/section assumptions, domain mapping, and covariate design.
Schema or expectation changes require a documented input-version change and a
new audit.

## Fit and result contracts

`run_screen` returns samples, exclusions, source/target DGELists, design, source
fit, target baseline/joint fits, mediator centering/scaling, complete results,
and the run signature. A failed mediator gate does not prevent fitting.

Result prefixes `a`, `c`, `cprime`, and `b` identify mediator diagnosis, baseline
target diagnosis, adjusted target diagnosis, and conditional M–Y effects.
Each has beta, moderated SE, CI, t, p, and gene-wide BH q. The a coefficient
is on the original source logCPM scale; b is target logCPM per one SD of M.
The saved mediator scaling is required to convert these scales. No a*b indirect
effect or indirect-effect p-value is estimated. The a-path repeats
across outcome rows; use `mediator_gates.tsv` for one row per screen.

`b_q` adjusts over all target genes in that fit. `b_q_global` adjusts across
all five primary M–Y families and is NA for subset runs. Sensitivity analyses
retain per-fit BH adjustment but are not included in primary pooled FDR.
`screen_hit`, `coefficient_shrinkage`, `higher_priority`, `direction_reversal`,
`attenuated_still_significant`, and `failed_gates` are separate fields.

Significance thresholds are strict < .05; loss is >= .05. Zero adjusted effects
count as nonreversing shrinkage. Relative shrinkage is NA for |c| <= 1e-8.
Signs conflicting with the original biological nomination are flagged, without
adding an additional gate. Identical source/outcome gene IDs in different
contexts are permitted and flagged.

## Provenance and cache rules

A SHA-256 signature includes configuration, all five screens, source-file
checksums/metadata, R source checksums, and package versions. Each cached voom
fit additionally hashes counts, normalization, design, donor block, and fitting
mode. Cache files are written through a same-directory temporary file and rename.
Changing any model-defining input prevents stale model reuse. The code commit
is recorded along with file fingerprints because working code may be uncommitted
and the remote repository may have unrelated changes.

Stages check primary artifact signatures before consuming them. Reports reject
stale result tables and mark older sensitivity/overlap statuses as unverified. After changing
R code, rerun the primary stage before sensitivities/reporting. Use separate
output directories for experiments or partial runs. Raw spot provenance is
recorded separately in the overlap stage because the large raw object is not
an input to the primary models.

## Operational failure handling

The audit stops on changed historical counts or matched sample expectations.
Primary fits record complete/failed status and abort the stage if any requested
screen fails. Unstable donor correlation and nonestimable coefficients are errors.
Sensitivity failures are recorded and make the job fail. The supplied launcher
runs sensitivity and overlap jobs independently after screening; the report job
requires both jobs to succeed. An overlap refit can return an explicit
`untestable_or_failed` status without terminating the entire overlap stage, so
the report and status table remain necessary for interpretation.

Shared-spot reruns may become untestable because of sample loss or mediator
filtering; these outcomes are recorded with reasons and are not counted as zero
hits. Review `overlap/status.tsv` and retention tables. Historical reproduction
reports numerical differences without replacing the historical results.
