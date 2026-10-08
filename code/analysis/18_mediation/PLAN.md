# Implementation plan and completion record

This document records the implemented scope and acceptance criteria for the
October 1, 2026 analysis. It is a consolidated plan record prepared with the
completed implementation; it is not a timestamped preregistration. Scientific
choices are specified in `config.R` and `screens.tsv`; rationale and limitations
are in [DESIGN.md](DESIGN.md) and
[mediation-analysis-decisions-notes.md](mediation-analysis-decisions-notes.md).

## Scope

Adapt the expression-screening framework from `LFF_spatial_ERC/code/22_Mediation`
to diagnosis-associated expression in spatialDLPFC_SCZ. The reference source
revision is `9d8d75d9df7e3052fd65ca558584ff7123cb3dfa`.
The new analysis occupies `code/analysis/18_mediation`, alongside the numbered
analysis folders. The eQTL/colocalization implementation is a separate workflow.

Five fixed hypotheses are evaluated:

| Source | Mediator | Target |
|---|---|---|
| Vascular | PTN | Neuropil |
| Vascular | PTN | Neuronal |
| Neuropil | FGF1 | Vascular |
| Neuropil | FGF2 | Vascular |
| Neuronal | FGF1 | Vascular |

Neuronal FGF2, genotypes, SNP PCs, expression PCs, and reference-project donor
exclusions are outside this implementation. No formal indirect-effect estimator
or causal signaling claim is part of the analysis.

## Phases and acceptance criteria

| Phase | Implementation and acceptance criteria | Completed evidence |
|---|---|---|
| Reference adaptation | Identify the screening rule, model outputs, variance treatment, and exact-baseline requirements; remove reference-project phenotype assumptions | `DESIGN.md`, `MODEL_SPECIFICATION.md`, `config.R` |
| Input inventory | Locate raw/count/logcounts objects, historical DE tables, metadata, labels, and stable gene identifiers | `reports/2026-10-01/audit/inventory.tsv`, input fingerprints |
| Sample alignment | Match donor–SpD keys, verify paired metadata, minimum spots, complete covariates, donor counts, and model rank | Audit sample/exclusion tables and `matched_designs.tsv` |
| Historical reference | Reproduce the original context-specific models and report numerical differences | `historical/reconciliation.tsv`; all three contexts match |
| Primary models | Fit all five source, exact target baseline, and mediator-adjusted models; retain full expressed-gene families and all gates | Five completed primary statuses; 40,329 pairs; four distinct baseline keys |
| Multiple testing and classification | Per-fit BH; pooled five-screen b BH; separate significance-loss, shrinkage, and direction flags | `primary/all_pairs.tsv.gz`, `mediator_gates.tsv`, independent verifier |
| Required sensitivities | Slide/RIN, matched logcounts, fixed weights for all screens; SpD07 exclusion for two neuronal screens | All 17 rows in `sensitivity/status.tsv` complete |
| Raw-spot audit | Reconstruct original GM pseudobulks exactly, quantify shared labels within matched sets | Raw count/spot reconciliation and overlap tables |
| Hit-triggered diagnostics | Whole-donor omission, between/within M, joint neuropil FGF1/FGF2, shared-spot removal; distinguish sample loss from zero hits | Implemented and synthetic-tested; no real-data triggers |
| Reporting and validation | Summaries, complete tables, five-page primary PDF, provenance, statistical and integration tests | `VALIDATION.md`, report archive, validation log |
| Review packaging | Impersonal documentation, explicit formulas and limits, full work log, portable archive, directory-only commit | This plan, `WORK_LOG.md`, archive manifest, review checks |

Paths in the evidence column below the first two rows refer to the archived run
under `reports/2026-10-01/` unless otherwise specified.

## Input availability and checks

Required inputs are defined as project-relative paths in `config.R`. The
executed JHPCE repository supplied all three pseudobulk objects, three historical
DE tables, donor metadata, and the raw kept-spot object. Existing inputs were
read-only. Local package availability was insufficient for the Bioconductor
analysis; execution used JHPCE `conda_R/4.5`.

Vascular–neuropil matching produced 272 observations from 61 donors (31 NTC,
30 SCZ). Vascular–neuronal matching produced 220 observations from 58 donors
(29 NTC, 29 SCZ). Primary GM domains are SpD02, SpD03, SpD05, SpD06, SpD07.
The data audit found no missing paired RIN values. Primary and slide/RIN baseline
matrices have ranks 8 and 24. Minimum sample requirements are enforced in code.

Historical vascular DE uses all seven SpDs; the matched screening models use
five GM domains. Historical and matched baseline eligibility therefore remain
separate. Reference-project baseline reuse is restricted to identical count,
normalization, sample, design, and donor-block inputs in this implementation.

## Inference and limitations

The primary model is `~ 0 + Dx + age + sex + SpD`, with donor blocking, adaptive
voom sample weights, and SCZ–NTC contrast. The adjusted target model adds M.
[MODEL_SPECIFICATION.md](MODEL_SPECIFICATION.md) defines the expression scale,
contrasts, moderation, sensitivity formulas, and coefficient units.

Nominal p < 0.05 is retained for historical and matched diagnosis gates. All four
historical mediator nominations fail gene-wide FDR < 0.05. Conditional M–Y tests
use full-family BH correction. Same-cohort selection, dependence across tests,
significance loss as a selection rule, mixed-cell measurements, shared spots,
and 58–61 independent donors limit interpretation. The composite rule has no
calibrated mediation-discovery FDR. A p-value crossing a threshold does not test
an indirect effect or a coefficient difference.

## Execution and failure behavior

The supplied launcher submits screening first, then independent sensitivity and
64 GB raw-spot jobs, then a report job dependent on both. The first completed run
used a combined primary/sensitivity/report allocation and a separate raw audit.
The work log records the execution difference and resource use.

Failed primary or required sensitivity fits cause stage failure. Shared-spot
refits record untestable/failed status and reasons. Zero hits from successful
fits remain distinct from failures. Complete stage statuses, run signatures,
and archived results are checked before interpreting the report.

## Completed outcome and remaining extensions

All five primary screens and all 17 required sensitivities returned zero
qualifying hits. All primary source mediator diagnosis gates failed. Conditional
M–Y associations are retained for inspection even where other gates fail.
Hit-triggered real-data diagnostics were not required under the fixed plan.

The following are outside the completed plan: independent replication,
functional experiments, formal indirect-effect estimation, an ordinary-limma
versus limma-trend comparison on an identical expression matrix, and a full
causal adjustment-set analysis. Their absence limits the claims supported by
this screen; it is not a pending failure of the implemented workflow.
