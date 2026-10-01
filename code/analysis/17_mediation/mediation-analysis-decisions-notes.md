# Analysis decisions and methodological basis

This record describes the implemented October 1, 2026 analysis. It distinguishes
fixed analysis choices, their statistical implications, completed diagnostics,
and extensions that have not been implemented. The full model specification is
in [MODEL_SPECIFICATION.md](MODEL_SPECIFICATION.md); execution evidence is in
[VALIDATION.md](VALIDATION.md).

## Decision record

| Dimension | Implemented choice | Purpose and limitation |
|---|---|---|
| Screening rule | Historical and matched diagnosis gates, conditional mediator–outcome FDR, and loss of adjusted diagnosis significance | Retains the significance-loss screen while reporting effect size and sign separately; it is not an indirect-effect test. |
| Observation unit | Matched donor–SpD pairs with donor blocking | Compares contexts within the same donor/domain; repeated domains are not independent donors. |
| Scope | Five fixed PTN/FGF screens; neuronal FGF2 excluded | Defines the testing family before fitting; exclusion does not establish biological absence. |
| Primary engine | TMM-normalized count input to adaptive-span voomLmFit with sample weights | Models observation-level precision; transformed logCPM is the fitted expression scale. |
| Primary covariates | Diagnosis, age, sex, SpD | Preserves the manuscript covariate set for the primary comparison. Sufficiency for causal adjustment is not established. |
| Expanded covariates | Slide and RIN as required sensitivity | Measures dependence on this expanded adjustment set; full rank establishes estimability only. |
| Outcome eligibility | Historical nominal DEG and matched-baseline nominal diagnosis association | Restricts screening eligibility while retaining full-gene fits and FDR families. Historical selection is not independent replication. |
| Failed mediator gate | Fit all five screens and retain the failed flag | Preserves complete association results without treating every association as a screening hit. |
| Multiplicity | Nominal diagnosis gates; BH FDR < 0.05 for conditional mediator–outcome association | BH applies to the stated association family, not the combined multistage screening rule. |
| Overlap | Audit both context pairs; remove shared spots and refit only for screens with primary hits | Assesses reuse of measurements; reaggregation also changes composition, sample availability, and precision. |

The count-based primary engine and manuscript covariate set are separate
choices. The historical reference analysis uses stored logcounts; it does not
become the primary analysis by sharing the covariate set.

## Expression scale and variance treatment

| Workflow | Fitted expression | Variance treatment |
|---|---|---|
| Ordinary limma on logCPM | Stored normalized logCPM | Gene-specific residual variances followed by default empirical-Bayes moderation; no abundance trend in the prior variance. |
| Limma-trend | Normalized logCPM | `eBayes(..., trend = TRUE)` models an abundance-dependent prior variance across genes. |
| Voom | Count-derived logCPM | Observation-level precision weights derived from a mean–variance trend; the implemented variant also estimates sample weights and accounts for donor blocking. |

The historical SPG workflow inspected for this adaptation stores TMM-normalized
log2 CPM in `logcounts`, then uses `lmFit()` with donor blocking and default
`eBayes()`. It supplies neither voom weights nor `trend = TRUE`. This description
applies to the audited SPG workflow and does not characterize other repository
analyses. Historical reproduction recovered all three nominal DEG counts and
numerically matched coefficients and p-values.

These workflows differ in variance treatment, not simply in whether expression
is log-transformed. Ordinary limma retains gene-specific residual variance
estimates. Donor blocking models repeated-observation correlation; it does not
replace mean–variance modeling. The original voom study compares these methods,
and the limma guide describes their operating conditions. The guide's library-
size guidance is not a measured property of this dataset or a universal cutoff.

Primary and matched-logcounts screens use the same matched sample sets and
primary-retained genes. They can differ in transformation, normalization history,
weights, correlation estimation, and residual degrees of freedom. Their
comparison therefore does not isolate weighting alone. An ordinary-limma versus
limma-trend comparison on an identical matrix, and systematic library-size and
mean–variance diagnostic summaries, are unimplemented extensions.

## Statistical interpretation

A transition from baseline p < 0.05 to adjusted p >= 0.05 does not test the
coefficient difference. It can reflect changed uncertainty or collinearity.
Coefficient shrinkage, sign reversal, and attenuation with persistent diagnosis
significance are separate output fields. No product-of-coefficients indirect
effect, indirect-effect confidence interval, causal proportion, or formal causal
mediation test is estimated.

BH adjustment uses all expressed outcome genes retained in each fit. The pooled
primary correction includes all 40,329 screen–gene pairs. Applying both corrections
does not calibrate the FDR of a list also selected through historical nomination,
nominal diagnosis gates, and significance loss. Selection and testing use the
same cohort. A failed nominal mediator gate describes this operational screen;
it does not establish absence of a biological pathway.

The independent biological replication unit is the donor. A consensus donor
correlation is used across repeated domains. This does not estimate every gene's
spatial covariance or separate between-donor from within-donor mediator effects.
The latter decomposition and whole-donor influence checks are implemented as
hit-triggered diagnostics, but were not triggered in the real-data run.

Slide/RIN sensitivity assesses a specified adjustment change. Additional
covariates do not automatically identify a causal adjustment set. Technical
quality, diagnosis, tissue composition, treatment, and other variables can have
relationships not resolved by these cross-sectional observations.

The raw-spot audit quantified nonexclusive SPG labels. Within matched samples,
28.00% of vascular spots were shared with neuropil and 20.85% with neuronal
labels. Removing these spots was not triggered by the zero-hit primary result;
its implementation was tested with synthetic data. No independent replication
or functional signaling experiment is part of this analysis.

## Completed results and scope of evidence

All five primary screens and all 17 required sensitivity analyses completed
with zero qualifying hits. Each primary source mediator diagnosis test failed
p < 0.05. The matched-logcounts analyses also failed this gate, so the difference
from historical nominal significance is not attributable solely to voom weights.
Full results remain available, including conditional expression associations
that fail other screening gates. The [validation record](VALIDATION.md) gives
sample counts, effect-test summaries, overlap, and completed checks.

## Methodological sources

References are identified by title and persistent identifier. These sources
support statistical principles; they do not establish a PTN/FGF mechanism in
this cohort. The software release documentation is mutable; the executed
package versions are recorded with the archived run.

1. [voom: precision weights unlock linear model analysis tools for RNA-seq read counts (2014)](https://doi.org/10.1186/gb-2014-15-2-r29).
2. [limma User's Guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf), RNA-seq sections on limma-trend and voom.
3. [The difference between “significant” and “not significant” is not itself statistically significant (2006)](https://doi.org/10.1198/000313006X152649).
4. [A comparison of methods to test mediation and other intervening variable effects (2002)](https://doi.org/10.1037/1082-989X.7.1.83).
5. [Controlling the false discovery rate: a practical and powerful approach to multiple testing (1995)](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).
6. [Circular analysis in systems neuroscience: the dangers of double dipping (2009)](https://doi.org/10.1038/nn.2303).
7. [Overadjustment bias and unnecessary adjustment in epidemiologic studies (2009)](https://doi.org/10.1097/EDE.0b013e3181a819a1).
8. [A practical solution to pseudoreplication bias in single-cell studies (2021)](https://doi.org/10.1038/s41467-021-21038-1). This supports a hierarchical-sampling principle rather than validation of the specific spatial model.
9. [A general approach to causal mediation analysis (2010)](https://doi.org/10.1037/a0020761).
10. [Bias in cross-sectional analyses of longitudinal mediation (2007)](https://doi.org/10.1037/1082-989X.12.1.23).
11. [limma reference manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf), entries for `eBayes`, `lmFit`, `duplicateCorrelation`, and `voom`.
12. [edgeR reference manual](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf), entries for `filterByExpr`, `calcNormFactors`, and `voomLmFit`.
