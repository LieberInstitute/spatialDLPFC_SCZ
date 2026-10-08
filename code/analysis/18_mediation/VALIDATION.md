# Execution and validation record

Date: 2026-10-01.

Result paths below refer to the [archived run](reports/2026-10-01/README.md).

**Complete:** all five primary screens, all 17 required robustness analyses,
historical reproduction, raw-spot auditing, reporting, and independent output
validation finished successfully. Every primary and robustness screen has zero
qualifying hits under the prespecified rule.

## Primary result

All five screens completed: **40,329 mediator–outcome pairs, zero qualifying
primary screening hits**. All five matched mediator diagnosis tests failed the
prespecified nominal p < 0.05 gate. This result applies to the configured screening
rule and is not evidence that the ligands have no biological role.

| Source mediator → target | Donors | Matched observations | Mediator Dx p | Source gene-wide BH FDR |
|---|---:|---:|---:|---:|
| Vascular PTN → neuropil | 61 | 272 | 0.574632 | 0.856660 |
| Vascular PTN → neuronal | 58 | 220 | 0.992975 | 0.996712 |
| Neuropil FGF1 → vascular | 61 | 272 | 0.257877 | 0.706140 |
| Neuropil FGF2 → vascular | 61 | 272 | 0.384778 | 0.796104 |
| Neuronal FGF1 → vascular | 58 | 220 | 0.343272 | 0.762488 |

The historical nominations used different sample/model definitions, including
all seven SpDs for vascular expression. Current primary screens use matched GM
samples and voom weighting. The matched-logcounts results below also fail the
source diagnosis gate, so the change from historical significance cannot be
attributed solely to voom weighting.

There are conditional mediator–outcome associations despite the failed source
diagnosis gates. Among genes satisfying both historical and matched-baseline
diagnosis eligibility, 385, 333, 70, 141, and 46 genes, respectively in the table
order above, have conditional mediator–outcome BH FDR < 0.05. These counts are
association results, not qualifying mediation-screen hits. All tested pairs
remain in `primary/all_pairs.tsv.gz` for inspection.

## Robustness results

All 17 required analyses completed without fit failures and returned zero
qualifying hits. The table gives source mediator diagnosis p-values; the full
gene-level results and FDR values are retained in `sensitivity/`.

| Source mediator → target | Slide + RIN | Matched logcounts | Excluding SpD07 |
|---|---:|---:|---:|
| Vascular PTN → neuropil | 0.339711 | 0.077511 | Not prespecified |
| Vascular PTN → neuronal | 0.687546 | 0.405449 | 0.949806 |
| Neuropil FGF1 → vascular | 0.234928 | 0.176222 | Not prespecified |
| Neuropil FGF2 → vascular | 0.346452 | 0.154333 | Not prespecified |
| Neuronal FGF1 → vascular | 0.203609 | 0.249517 | 0.431996 |

The five fixed-weight comparisons retain the primary source mediator fit, so
their source diagnosis p-values equal the primary values by construction; their
target comparisons also yielded zero hits. After SpD07 exclusion, both neuronal
screens retained 210 observations from 58 donors. Re-filtering within the primary
gene universes retained 12,509 neuronal outcome genes for PTN and 5,097 vascular
outcome genes for FGF1.

The closest matched-logcounts source result is PTN → neuropil (p = 0.0775).
These sensitivities assess specified model choices within the same cohort; they
provide neither independent replication nor evidence that the ligands lack a
biological role.

## Input and historical-model checks

- Matched donor–SpD counts and all candidate gene identifiers agree with the
  pre-implementation audit.
- Primary and slide/RIN model matrices have ranks 8 and 24, respectively.
- Historical nominal DEG counts reproduced exactly: 440 vascular, 1,789 neuronal,
  and 1,669 neuropil.
- Maximum absolute historical coefficient differences were below 9e-14;
  maximum p-value differences were below 2.1e-12.
- The raw spot object reconstructed the original GM pseudobulk counts exactly
  for every retained input gene: vascular 5,215 genes/272 samples, neuronal
  16,768/259, and neuropil 11,834/315. Contributing-spot counts also matched.

## Shared-spot audit

Fractions below use each primary screen's matched sample keys, not all raw spots
in the study. The reverse FGF screens use the same corresponding matched sets.

| Matched contexts | Shared spots | Vascular spots | Other-context spots | Shared/vascular | Shared/other |
|---|---:|---:|---:|---:|---:|
| Vascular–neuropil | 5,241 | 18,721 | 122,449 | 28.00% | 4.28% |
| Vascular–neuronal | 3,302 | 15,836 | 47,089 | 20.85% | 7.01% |

The audited SPG labels are nonexclusive. Real-data overlap-removal
refits, donor-influence diagnostics, between/within-donor diagnostics, and the
joint FGF model were not triggered because there were no primary hits. Their
implementation was exercised separately with synthetic data.

## Automated validation

Passed checks include:

- Shuffled sample matching; inconsistent diagnosis and duplicate-key rejection.
- Missing-covariate, constant-mediator, and rank-deficient-design rejection.
- SCZ–NTC contrast direction, a known synthetic mediator signal, and an
  independent negative-control mediator.
- Correct separation of significance loss, coefficient shrinkage, and sign reversal.
- Synthetic hit-triggered whole-donor influence and between/within-donor models.
- Synthetic joint FGF modeling and full-gene testing families.
- Exact synthetic raw-count reconstruction, removal of shared spots from both
  contexts, explicit insufficient-spot status, raw-count mismatch rejection,
  and zero-hit control flow. Forced triggers remain confined to temporary
  synthetic fixtures and are not scientific findings.
- Scheduler dependency and spaced-path dry-run checks, without submitting jobs.
- Independent Python verification of all 40,329 real-data rows: unique IDs,
  expected samples/genes, identical shared baselines, and per-screen and pooled
  BH corrections recomputed from exported p-values.
- Independent verification of all 17 robustness tables: complete tasks, matching
  run signatures, unchanged or correctly restricted gene universes, sample
  counts, and full-gene BH corrections recalculated from exported p-values.
- All five primary diagnostic PDF pages rendered and visually checked for
  readable labels, correct screen/sample annotations, and complete panels.

## Provenance and execution

R model signature:
`d7fb6893ce5cf7afaa6a2d627eb7359ca2bfb80766e0ec811c96cb589953be8b`.

Runtime: JHPCE `conda_R/4.5`, R 4.5.0 patched, edgeR 4.6.2, limma 3.64.0.
Full input/code fingerprints and session information are retained in the run's
`provenance` directory. Scientific R code was unchanged during this run.

- Full-data primary/robustness/report run: Slurm job `36114706`, completed with
  exit code 0 in 2 h 2 min 8 s; approximately 9.6 GiB peak resident memory within
  its 16 GB allocation.
- Raw-spot audit: job `36115731`, completed successfully; approximately 24 GiB
  peak resident memory within its 64 GB allocation.
- Final synthetic integration fixture: job `36115761`, completed successfully.
- Primary report preview for visual QA: job `36115782`, completed successfully;
  the full-data job subsequently regenerated the final report with all 17
  robustness tasks marked complete.

This first execution used a combined primary/robustness job. The supplied
launcher now separates screening, sensitivities, and spot auditing, allowing
independent checks to run concurrently and producing the final report only
after both robustness jobs succeed.
