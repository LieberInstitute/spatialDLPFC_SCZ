# Design decisions and interpretation

## Biological scope

The analysis evaluates diagnosis → vascular PTN → neuropil or neuronal
expression, and diagnosis → neuropil FGF1/FGF2 or neuronal FGF1 → vascular
expression. Historical nominal SCZ DEGs define outcome eligibility. The five
screens are fixed in `screens.tsv`.

This is independent of the eQTL/coloc pipeline. Genotypes, SNP PCs, expression
PCs, APOE carrier contrasts, and ERC donor exclusions are not inputs.

## Data and matching

Use the three original count/logcounts SpatialExperiment objects in
`processed-data/rds/PB_dx_spg`, not donor-only eQTL exports. Keys are donor and
SpD, with section identity checked. Every current donor contributes one section.
Multiple sections per donor are rejected because no aggregation rule is specified.
The primary GM domains are spd02, spd03, spd05, spd06, and spd07.

The vascular object contains 5,215 genes/352 columns/61 donors; neuronal contains
16,768/308/63; neuropil contains 11,834/427/63. Vascular–neuropil matching yields
272 observations and 61 donors, whereas vascular–neuronal matching yields
220 observations and 58 donors. Primary and slide/RIN designs have ranks 8 and
24, respectively. There are no missing paired RIN values in the audited inputs.

Use `test_SPD_pseudo_vasc_pos.csv`, `neun-dx_DEG-GM.csv`, and
`neuropil-dx_DEG-GM.csv` as the historical selection sources. Their nominal gene
counts are 440, 1,789, and 1,669. Historical vascular models use all seven SpDs;
new vascular outcome models use the matched five-domain sample sets. These analyses have different sample sets and are reconciled separately.

## Model design

Source and target counts are filtered separately with `filterByExpr` using the
matched baseline design, then normalized with TMM after resetting library sizes.
Source mediator values are taken from the source baseline voom EList and scaled
by the observed matched-sample SD. The raw scale and centering are saved.

Baseline and source models: `~ 0 + Dx + age + sex + SpD`.
Joint target model: `~ 0 + Dx + age + sex + SpD + M`.
Diagnosis contrast: SCZ minus NTC. Donor is the repeated-measures block.
Each fit uses adaptive-span voom with sample weights and limma empirical Bayes.
Source diagnosis testing is performed in the full source gene universe.

Keep target genes, normalization, and samples identical within each nested
comparison. Only identical model inputs permit baseline reuse. Four target
baselines suffice for the five primary screens; neuropil FGF1/2 share one.
A baseline from a different sample set is not reused.

Historical DEG eligibility and matched baseline nominal significance are both
required. All genes are nevertheless fitted and included in each M–Y FDR family.
All five source candidates are tested even if the nominal a-path gate fails.
The expected primary target universes total 40,329 M–Y tests; actual retained
universes and denominators are reported.

## Robustness

- Add slide and RIN using the same samples and genes.
- Reproduce manuscript models on stored logcounts and original sample sets;
  the historical correlation design excludes SpD, matching upstream code.
- Compare matched-sample logcounts results using the coherent full covariate design.
- Hold baseline voom expression/weights and a newly estimated baseline donor
  correlation fixed in both nested models, isolating inclusion of the mediator
  from re-estimation of voom weights/correlation. This diagnostic uses lmFit;
  it is not an exact recreation of voomLmFit's sparse-count degrees of freedom.
- For hit pairs, leave out whole donors and report coefficient changes without
  treating the selected-gene diagnostic refits as new inferential tests.
- Decompose M into its donor mean and within-donor deviation, fit the full gene
  universe, and label hit pairs. Means use available matched SpDs, so uneven
  domain coverage still limits this diagnostic.
- Exclude spd07 in neuronal-involving screens: its primary matched coverage is
  only seven NTC and three SCZ observations. Re-filter within each primary gene
  universe after removing these observations; report any genes lost to filtering.
- If neuropil FGF hits exist, include FGF1 and FGF2 together; report collinearity.
- Require the raw spots to reconstruct the existing GM pseudobulk counts and
  contributing-spot counts exactly, guarding against input-version differences.
  Audit overlap from the actual raw spot labels. For screens with hits, remove
  shared spots from both contexts, reaggregate counts, require ten remaining
  spots, re-match donors/SpDs, re-filter within the primary gene universes, and
  refit. Retain the original nominations and historical outcome lists.

## Caveats and scientific boundaries

The historical source nominations do not pass FDR < 0.05: historical
FDR values are .407 for vascular PTN, .330 for neuropil FGF1, .308 for neuropil
FGF2, and .460 for neuronal FGF1. Nominal selection does not
provide genome-wide evidence that diagnosis changes these mediators.

Selection and screening reuse the same cohort. M–Y BH corrections describe
those association families; they do not establish a mediation-discovery FDR
for the combined selection, gate, and attenuation procedure.

A p-value becoming nonsignificant is not a test of coefficient attenuation.
It can result from larger uncertainty or collinearity; a coefficient can also
shrink without losing significance. Report magnitude, direction, uncertainty,
and the original screening flag separately. Relative shrinkage is not a causal
proportion mediated, and coefficient differences are not assigned causal p-values.

Cross-sectional postmortem observations do not establish temporal ordering.
Clinical diagnosis is not randomized. Medication, smoking, illness duration,
cellular composition, tissue quality, and other common causes can influence
both mediator and outcome. Slide/RIN sensitivity cannot remove all confounding.
RIN adjustment is itself a modeling sensitivity, not proof of a causal adjustment set.

There are 58–61 independent donors, not 220–272 independent samples. Blocking
uses a consensus correlation and cannot model every spatial covariance pattern.
Missing context/domain measurements and variable spot counts limit generalization.
SPGs contain mixtures of cells and can share spots, inducing mechanical
expression correlation. Removing overlap also changes tissue composition and
precision; a failed rerun must be distinguished from loss of association.

The available pseudobulks have already undergone upstream expression filtering,
especially the vascular object. Missing genes are untested, not biologically absent.
The five univariate mediator screens cannot establish unique causal pathways
among correlated ligands. A robust result is a hypothesis for functional follow-up.

Methodological references:
- Cross-sectional mediation bias (2007): https://pubmed.ncbi.nlm.nih.gov/17402810/
- Causal identification assumptions: https://pmc.ncbi.nlm.nih.gov/articles/PMC3659198/
- edgeR documentation: https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf
