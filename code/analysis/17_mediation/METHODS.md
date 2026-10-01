# Draft methods: exploratory PTN/FGF expression screening

Five prespecified cross-microenvironment expression hypotheses were evaluated:
vascular PTN with neuropil and neuronal outcomes; neuropil FGF1 and FGF2 with
vascular outcomes; and neuronal FGF1 with vascular outcomes. Source and target
pseudobulks were matched by donor and gray-matter spatial domain (SpD02, SpD03,
SpD05, SpD06, and SpD07). Pseudobulks required at least ten contributing spots.
Neuronal FGF2 was not included.

Within each matched analysis, source and target counts were separately filtered
using edgeR filterByExpr with the baseline design and TMM-normalized. Models used
edgeR voomLmFit with an adaptive mean–variance trend, sample weights, donor blocking,
and limma empirical-Bayes moderation. Baseline models included diagnosis, age,
sex, and spatial domain. Diagnosis effects were expressed as SCZ minus NTC.
Source mediator expression was obtained from the corresponding normalized voom
logCPM matrix and standardized to one SD over matched observations. One mediator
was added to each target model. Baseline and adjusted comparisons used identical
target samples, genes, and normalization.

A pair qualified for exploratory screening when the target gene was a historical
nominal SCZ DEG (p < .05), retained a nominal diagnosis association in the matched
voom baseline (p < .05), the source mediator had a nominal matched diagnosis
association (p < .05), the conditional mediator–outcome association passed BH
FDR < .05 across all expressed target genes in that screen, and the adjusted
target diagnosis association was no longer nominally significant (p >= .05).
All expressed genes were fitted and included in the relevant testing families.
All five candidate screens were fitted regardless of the mediator gate result.

Coefficient shrinkage and direction changes were reported separately from loss
of statistical significance. Higher-priority screening hits additionally required
same-direction coefficient shrinkage and mediator–outcome BH FDR < .05 across
all five primary screens. This prioritization does not provide a calibrated
mediation-discovery FDR for the full selection procedure. Relative changes in the
diagnosis coefficient were not interpreted as causal proportions mediated.

Sensitivity analyses added slide and RIN, compared matched logcounts models,
held baseline weights/correlation fixed in nested models, and excluded sparse
neuronal SpD07 observations. Hit-triggered diagnostics were implemented for whole-donor
influence, between-/within-donor mediator components, joint neuropil FGF1/FGF2
modeling, and removal of shared spots from both contexts. Shared-spot counts were
audited for both context pairs. No real-data hit-triggered refits were required
because the primary screens returned zero hits. Synthetic fixtures exercised
these code paths, including reaggregation and expression filtering within the
primary gene universes.

These analyses are observational and cross-sectional. Historical mediator
nominations were nominal and did not pass gene-wide FDR correction. Selection
and testing reused the same cohort. Mixed-cell measurements, shared spots,
unmeasured confounding, uncertain temporal ordering, and limited independent
donor numbers preclude causal or functional claims. Significance loss alone is
not a test of coefficient change. Results identify patterns for follow-up rather
than demonstrating ligand-mediated signaling.

The executed run used 61 donors/272 matched observations for vascular–neuropil
comparisons and 58 donors/220 observations for vascular–neuronal comparisons.
All 40,329 primary pairs failed the combined screening rule; all five source
diagnosis gates failed (p = 0.258–0.993). All 17 required sensitivity analyses
completed with zero qualifying hits. Software versions were R 4.5.0 patched,
edgeR 4.6.2, and limma 3.64.0.
