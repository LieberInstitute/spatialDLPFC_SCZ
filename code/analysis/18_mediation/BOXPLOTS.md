# Candidate expression boxplots

The complete nominal screen is displayed: 110 mediator–outcome pairs, 90 outcome
genes and five screens. A biological shortlist contains 11 genes and all 16
qualifying pairs involving them. Selection changes neither the screening gates
nor the fitted models. All unselected pairs remain in the complete atlases.

## Plot construction

The ERC `plotMedComp` approach in
`LFF_spatial_ERC/code/22_Mediation/02_E2_E4_med_plots.Rmd` compares outcome
expression before and after mediator adjustment, followed by mediator expression
in its source context. ERC uses `jaffelab::cleaningY` to retain APOE group effects
while removing other design effects. The SCZ adaptation retains diagnosis and
uses the saved donor-blocked GLS coefficients rather than an unweighted refit.

For expression vector y, nuisance design Z and its saved coefficient vector b,
plotted values are `y - (Z - columnMeans(Z)) b`:

- Baseline outcome: remove age, sex and SpD contributions from the baseline fit.
- Adjusted outcome: remove age, sex, SpD and standardized mediator contributions
  from the joint fit.
- Source mediator: remove age, sex and SpD contributions from its source fit,
  using the same matched observations as the outcome panels.

Centering preserves the original mean expression and places both outcome panels
at the same average covariate profile. The two outcome panels share y-limits;
the mediator uses its own scale. Nuisance coefficients differ between baseline
and joint fits, so the middle panel is not simply the first panel minus M.
Values use the saved **logcounts** assay, not a newly computed logCPM transform.
These are adjusted observational distributions, not counterfactual outcomes.

Two complete versions are provided:

- **Donor–SpD:** each point is one matched observation. Boxes describe all
  observations; repeated observations from a donor are explicitly identified
  in captions and group counts.
- **Donor means:** adjusted values are averaged across each donor's available
  matched SpDs, yielding one point per donor and equal donor weight in the
  displayed distribution. Domain coverage can differ between donors. This is
  a descriptive companion, not a balanced-domain estimate or a new model.

Boxes show medians and interquartile ranges, with whiskers to the furthest values
within 1.5 IQR. Every observation is shown, including outliers; jitter is fixed
and aligned across panels. Annotated coefficients, nominal p-values and BH FDR
are copied from the unchanged full-family limma results. The donor-mean version
retains those original model statistics; no t-tests or other tests are performed
on the plotted values. The conditional mediator coefficient is per mediator SD.
Diagnosis attenuation is a coefficient change, not a proportion mediated.

The source panel shows the matched mediator association, including its
non-significant diagnosis p-value. The accepted screening premise uses the
original manuscript ligand nomination, not a matched Step-2 significance gate.
Loss of diagnosis significance after adjustment does not establish mediation.
Shared cell composition, overlapping spots, selection effects and confounding
remain possible explanations for the associations.

## Biological shortlist

The current [SCZ PNN manuscript](https://docs.google.com/document/d/1XZjI8Q4ArNRrdcBLz_eWI66ORVSbbuzTTClkq48l_lY/edit)
was read on 2026-10-21. Results 2.6 describes synaptic/ECM, vascular and metal-ion
homeostasis themes and PTN/FGF trophic axes. Results 2.7 and Figure S42 discuss
LY6H; Results 2.4 names GAP43. These manuscript connections motivate display
selection, not additional evidence that the nominated ligand mediates an outcome.

| Gene | Pairs | Selection rationale |
| --- | ---: | --- |
| LY6H | 1 | Strongest direct manuscript connection: neuronal SCZ expression and spatial eQTL/colocalization; PTN → neuronal. |
| TIMP2 | 1 | Matrix-metalloproteinase regulation connects to the ECM/PNN theme; PTN → neuropil. [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/7077). |
| IGSF9B | 1 | Synaptic adhesion and inhibitory-synapse biology; PTN → neuropil. [Developmental study](https://pubmed.ncbi.nlm.nih.gov/23751499/); [adult mouse localization](https://pmc.ncbi.nlm.nih.gov/articles/PMC12675716/) is not exclusively inhibitory. |
| GAP43 | 2 | Named in the manuscript for neuronal structure/outgrowth; compare FGF1 and FGF2 neuropil → vascular. |
| APOLD1 | 1 | Endothelial junction/permeability context; FGF1 neuronal → vascular. [Functional study](https://pubmed.ncbi.nlm.nih.gov/35638551/). |
| CNNM2 | 1 | External cortical/synaptic SCZ evidence; PTN → neuropil. [Human and experimental study](https://pubmed.ncbi.nlm.nih.gov/37715107/). |
| RGS4 | 1 | Prior cortical SCZ expression work, with [positive](https://pubmed.ncbi.nlm.nih.gov/11326297/) and [negative](https://pubmed.ncbi.nlm.nih.gov/18622782/) reports; FGF1 neuronal → vascular provides a context comparison. |
| ADCYAP1R1 | 1 | PACAP receptor signaling; FGF1 neuronal → vascular. It is not an FGF receptor. [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/117). |
| MT3 | 2 | Representative of the MT1G/MT1M/MT3-driven metal-response GO pattern; compare neuropil FGF1 and FGF2. |
| AGT | 2 | Shared FGF1/FGF2 vascular outcome contributing to the growth GO term; useful ligand comparison. |
| FAM107A | 3 | Present in all three FGF → vascular screens; compare ligand and source microenvironment. |

The shortlist spans all five screens and includes biological and cross-screen
comparisons rather than only the smallest p-values or largest attenuation.
`boxplot_shortlist.tsv` records the rationale and sources; `selected_pairs.tsv`
adds page numbers, exact pair membership and screening statistics. A vascular
SPG label does not establish endothelial cellular origin. These five screens
contain no PNN-target model, so PNN/ECM relevance is thematic. The shortlisted
outcomes are not asserted to be direct PTN or FGF receptors.

## Files and reproduction

Plots: `plots/18_mediation/nominal_p05/candidate_boxplots/`:

- `candidate_boxplots_all_110.pdf`: complete 110-page atlas.
- `candidate_boxplots_all_110_donor_means.pdf`: complete donor-mean companion.
- `candidate_boxplots_selected.pdf`: 16 selected pairs in biological-priority order.
- `candidate_boxplots_selected_donor_means.pdf`: the same ordered subset, donor means.
- `selected/`: individual vector PDF/SVG and 300-dpi PNG for each selected pair.

Data: `processed-data/18_mediation/run_2026-10-08_limma/candidate_boxplots/`:
`plot_values.tsv.gz`, `pair_plot_index.tsv`, `selected_pairs.tsv`,
`model_verification.tsv` and `provenance.json`. Plotting values contain only
anonymous donor codes, diagnosis, SpD and the three adjusted expression values;
original donor identifiers and demographic covariates are omitted. Donor codes
are consistent across screens but are not accompanied by an identifier map.

From the project root, using the existing figure requirements:

```bash
python code/analysis/18_mediation/plot_candidates.py \
  --source-run processed-data/18_mediation/run_2026-10-08_limma \
  --plots plots/18_mediation/nominal_p05/candidate_boxplots
python code/analysis/18_mediation/tests/verify_candidate_boxplots.py \
  --source-run processed-data/18_mediation/run_2026-10-08_limma \
  --plots plots/18_mediation/nominal_p05/candidate_boxplots
```

Verification uses Poppler `pdfinfo` and `pdftotext`. Regenerating the plotting
snapshot additionally requires the five original `primary/<screen>.rds` fitted
artifacts, retained outside Git:

```bash
Rscript code/analysis/18_mediation/export_candidate_boxplots.R \
  processed-data/18_mediation/run_2026-10-08_limma
```

No Bioconductor installation is needed for extraction: the saved fit components
are read as ordinary lists. Base-R reconstruction of the GLS coefficients using
the saved donor correlation agrees within 1.2e-12 for all plotted outcomes and
mediators across five screens. Checks also verify sample alignment, mediator
standardization, result coefficients, centered-value means and input membership.
PDF verification checks all 252 atlas pages plus 16 individual PDFs against the
pair indexes, and verifies input/output hashes. Automated layout bounds checks
cover every generated panel; representative observation and donor-mean pages
were visually inspected. These checks verify plotting fidelity, not causal validity.
