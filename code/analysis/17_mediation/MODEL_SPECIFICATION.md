# Model specification and inference

## Observation and variable definitions

An observation is a donor–SpD pseudobulk available in both source and target
contexts. Matching uses the composite donor/domain key and verifies section,
diagnosis, age, sex, slide, and RIN agreement. The inputs contain one section per
donor. Multiple sections require a revised aggregation rule; they are rejected
by the current loader. At least ten contributing spots are required in each
context. Complete slide/RIN metadata are required even for the primary analysis,
so the required expanded sensitivity uses the same observations.

`Dx` is a factor with levels NTC and SCZ. `sex` and `SpD` are factors; `age`
is numeric. The expanded design treats `slide_id` as a factor and `rin` as
numeric. Age and RIN are not centered or standardized by the implementation.
Factors use R's model-matrix contrasts; actual column names define the tested
contrast. Every design must be full rank, have at least five residual design
degrees of freedom, at least ten donors total, and at least five donors per
diagnosis. The current matched sets exceed these minima.

The archived [design column table](reports/2026-10-01/report/design_columns.tsv)
records the column names, order, sample counts, and ranks reconstructed from
archived primary samples using `make_design()`.

## Expression preparation

Source and target count matrices are prepared separately on exact matched
samples. `filterByExpr()` uses the baseline design. After filtering, library
sizes are reset to retained count sums and TMM normalization factors are
recomputed. Input pseudobulks already have upstream gene filtering; the current
analysis cannot recover discarded genes.

The source baseline fit supplies the source mediator's normalized logCPM vector
`L`. Its observed matched-sample mean and sample SD define
`M = (L - mean(L)) / sd(L)`. Scaling is over donor–SpD observations rather than
equally weighted donor means. A missing, filtered-out, or constant mediator
causes a fit failure. Saved fit artifacts include the scaling parameters.

## Primary formulas

Source diagnosis model and target baseline model:

```r
~ 0 + Dx + age + sex + SpD
```

Target mediator-adjusted model:

```r
~ 0 + Dx + age + sex + SpD + M
```

For donor d and domain s, the baseline mean can be written as:

```text
E[Y_ds] = beta_NTC I(Dx_d = NTC) + beta_SCZ I(Dx_d = SCZ)
          + beta_age age_d + sex terms + SpD terms
```

The adjusted model adds `b M_ds`. The source model has the same baseline
structure with source expression as its response. `0 + Dx` provides separate
diagnosis columns, so the diagnosis contrast is always `DxSCZ - DxNTC`.
The primary baseline matrix has rank 8; adding M adds one column when estimable.

Donor enters the fitting call as a repeated-observation block, not as a column
of donor fixed effects:

```r
edgeR::voomLmFit(
  dge, design = design, block = samples$donor,
  adaptive.span = TRUE, sample.weights = TRUE, keep.EList = TRUE
)
```

This call models the mean–variance relationship, observation precision, sample
weights, and within-donor correlation using the installed edgeR implementation.
There is one consensus within-donor correlation in each fit. The implementation
does not fit a separate donor correlation for every gene or an arbitrary spatial
covariance matrix. Source, baseline, and adjusted fits estimate their variance
components separately. Baseline and adjusted targets use identical count data,
retained genes, TMM normalization, and observations.

## Contrasts, moderation, and units

Each coefficient table is obtained with `contrasts.fit()`, followed by default
`eBayes()` and `topTable(..., number = Inf, sort.by = "none", confint = TRUE)`.
No `trend = TRUE` or `robust = TRUE` argument is supplied to this moderation
step. Primary variance modeling comes from the voom fit. The exported standard
error uses the unscaled coefficient standard error times the square root of the
moderated posterior variance; confidence limits and p-values use moderated
inference.

| Prefix | Test | Units |
|---|---|---|
| `a` | Source mediator SCZ–NTC contrast | Source logCPM difference; not standardized M units |
| `c` | Target baseline SCZ–NTC contrast | Target logCPM difference |
| `cprime` | Target adjusted SCZ–NTC contrast | Target logCPM difference conditional on M |
| `b` | M coefficient in adjusted target model | Target logCPM per one observed SD of source mediator |

The a-path table is fit across all retained source genes before extracting the
mediator. Its p-value is repeated across outcome rows. Its gene-wide q-value
therefore refers to the source universe, not the number of repeated rows.
`a_q_five_candidates` additionally adjusts the five nominated source tests.
No `a*b` estimate is computed; the a and b output scales differ.

## Testing families and screening flags

Per-fit BH q-values are calculated separately for c, cprime, and b across the
retained outcome universe. Only b q-values are used in the mediator–outcome gate.
The primary pooled b correction covers all five screen–gene families, totaling
40,329 tests. A gene appearing in multiple screens contributes one test per
screen. Subset runs have no complete pooled primary correction. Sensitivity
families are adjusted separately and are excluded from the primary pooled family.

A screening hit is the conjunction of:

1. Historical target diagnosis p < 0.05.
2. Matched target baseline diagnosis p < 0.05.
3. Matched source mediator diagnosis p < 0.05.
4. Conditional mediator–outcome BH q < 0.05 within the screen.
5. Adjusted target diagnosis p >= 0.05.

The additional `higher_priority` label requires `abs(cprime) < abs(c)`, no sign
reversal, and pooled b q < 0.05. Zero adjusted effect counts as nonreversing
shrinkage. Relative shrinkage is `1 - cprime/c`, omitted for `abs(c) <= 1e-8`.
This is a descriptive coefficient comparison, not a causal proportion. Expected
mediator direction and identical source/target gene IDs are flags, not extra
eligibility gates. All fitted rows remain available regardless of gate status.

## Reference and sensitivity formulas

| Analysis | Design and variance treatment | Sample/gene family |
|---|---|---|
| Historical reference | Stored logcounts; `duplicateCorrelation()` with `~ 0 + Dx + age + sex`; `lmFit()` with `~ 0 + Dx + age + sex + SpD`; default moderation | Original seven-domain vascular or five-domain neuronal/neuropil sets; reproduces upstream model behavior |
| Matched logcounts | Stored logcounts; baseline and adjusted full designs as above; a new donor correlation estimated for each fit; default moderation | Primary samples and genes; M restandardized from stored source logcounts |
| Slide/RIN | Add `slide_id + rin` to source/baseline and adjusted designs | Primary samples and genes; primary standardized M retained; baseline rank 24 |
| Fixed weights | Reuse target baseline EList expression and weights; estimate one baseline donor correlation and hold it fixed in `lmFit()` baseline/adjusted fits | Primary samples and genes; primary source fit and M retained |
| Without SpD07 | Primary voom formulas after removing SpD07; re-filter and renormalize within primary source and target gene universes | Two neuronal-involving screens; 210 observations/58 donors; M restandardized after refitting |

The fixed-weight diagnostic uses ordinary `lmFit()` residual degrees of freedom;
it does not reproduce every sparse-count degrees-of-freedom adjustment in
`voomLmFit()`. Changes between the historical and primary models include sample
matching, gene filtering, normalization, correlation design, and weighting.
Neither comparison alone attributes all differences to a single component.

## Hit-triggered diagnostics

- Whole-donor omission removes every domain for that donor. Target baseline
  expression, weights, and baseline donor correlation are retained; coefficients
  are exported without treating these selected-gene refits as independent tests.
- Between/within decomposition uses the donor mean of M and `M - donor mean`.
  Both terms enter a full-gene target model using fixed baseline weights and
  correlation. Unequal domain coverage can affect donor means.
- Joint neuropil FGF1/FGF2 modeling includes both standardized mediators, verifies
  identical matched keys, and reports their correlation and design condition
  number. Each mediator coefficient uses its full-gene testing family.
- Shared-spot removal deletes spots labeled in both contexts from both count
  aggregations. At least ten remaining spots are required. Samples are rematched,
  genes re-filtered within the primary universes, and voom models refitted.
  Original nominations and historical eligibility are retained. Sample loss or
  failed estimation is recorded separately from a successful zero-hit fit.

None of these hit-triggered models ran on the real-data cohort because there
were no primary hits. Synthetic integration fixtures exercise their control
flow and failure handling. The raw overlap audit and count reconstruction ran
on the real data regardless of hit count.

## Review boundaries

The procedure estimates conditional expression associations under specified
models. Diagnosis is observational; time order and unmeasured common causes are
unresolved. Mixed-cell pseudobulks and shared spots can induce association.
Historical nomination and testing reuse the same cohort. The composite screening
rule has no established mediation-discovery FDR. Full model details and these
limitations apply equally to positive and zero-hit outputs.
