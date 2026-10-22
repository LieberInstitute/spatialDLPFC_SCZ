# Nominal mediation: overlap and GO enrichment

The unchanged October 8 fits yield 110 candidate mediator–outcome pairs,
90 unique Ensembl outcome genes, five screens, four mediator/source combinations,
and three ligand genes. These are exploratory screening results, not established
causal mediation. The foreground is identical in both GO analyses.

## Results

| Screen | Candidates | Eligible-DEG background | All-tested background |
| --- | ---: | ---: | ---: |
| PTN vascular → neuropil | 39 | 1,187 | 11,834 |
| PTN vascular → neuronal | 22 | 1,293 | 16,768 |
| FGF1 neuropil → vascular | 17 | 258 | 5,215 |
| FGF2 neuropil → vascular | 26 | 258 | 5,215 |
| FGF1 neuronal → vascular | 6 | 200 | 5,215 |

Background counts precede Entrez mapping and ontology-specific annotation filtering.

**Overlap:** PTN target sets share CAMK1G and ZNHIT1. The three FGF screens
contain 31 unique vascular outcomes; the two neuropil FGF screens share 16,
including FAM107A and SLC24A2 also found in the neuronal FGF1 screen.
Exact intersections total 90 genes and 110 pair incidences. Matching uses Ensembl
IDs; the same gene across different targets denotes gene identity, not the same
expression measurement.

**Eligible-DEG GO:** no term reaches within-family BH p < .05 across 3,143
query-overlapping term tests in 15 screen/ontology families. The minimum is
.06396 for plasma membrane and cell periphery in FGF1 neuronal → vascular.
This tests functional overrepresentation after conditioning on the two DEG gates.

**All-tested GO:** 11 screen–term associations (10 distinct terms) reach
within-family BH p < .05 across 4,560 tests:

- FGF2 neuropil → vascular: nine biological-process terms, including growth and
  metal-ion responses. Seven metal-related terms share exactly MT1G, MT1M and MT3;
  these are overlapping annotations, not seven independent biological findings.
- FGF1 neuropil → vascular: regulation of growth, six genes, BH p = .01719.
- PTN vascular → neuronal: protein serine/threonine kinase binding, NAIP and
  CCDC102B, BH p = .03101.
- No cellular-component term passes in any screen under either background.

Eight associations survive pooling the 15 query-overlap testing families.
Seven survive the separate sensitivity correction over all size-eligible terms,
including zero-overlap terms at p = 1; all seven are in FGF2 neuropil → vascular.
These are separate corrections, not a combined guarantee. The all-tested
analysis characterizes candidates relative to the measured target transcriptome;
it can include functional structure inherited from DEG selection.

## Approach and interpretation

The ERC reference uses a pooled FDR < .05 DEG universe from Oligo.3 and Astro.1–3.
The eligible-DEG analysis follows that restriction in principle, but uses the
actual outcome pool separately for each screen: manuscript diagnosis p < .05
and matched-baseline diagnosis p < .05. Source and target genes are not pooled.
The all-tested analysis is a complementary biological characterization. Both
backgrounds are retained; significance does not determine which is reported.
This distinction follows the [GO reference-list guidance](https://geneontology.org/docs/go-enrichment-analysis/)
that the reference should represent the population from which the analyzed list
was selected, with the conditioning stage stated explicitly.

- Candidate gates remain historical p < .05, matched baseline p < .05,
  mediator coefficient p < .05, and adjusted diagnosis p ≥ .05. All original
  92 pairs are retained, with 18 additions. FDR columns remain annotations.
- Matched mediator–diagnosis p-values range from .0679 to .3961 and are not
  selection gates. The original manuscript association supplies this screening
  premise. This is an explicit compromise, not a fully satisfied matched-data
  Baron–Kenny sequence.
- All 110 pairs have same-sign coefficient attenuation (5.82–37.13%), without
  adding an attenuation gate. Crossing the significance threshold does not test
  an indirect effect or prove complete mediation. No indirect-effect confidence
  intervals or causal identification are supplied here.
- Donor-blocked fitted models use 58 or 61 donors, with 220 or 272 observations;
  observations are not independent donors. GO analyzes gene sets and does not
  provide new donor-level evidence.
- ORA assumes exchangeability under its gene-sampling null. DEG conditioning
  does not remove expression/detection bias, gene correlation, selection effects,
  or possible cell-composition confounding. No empirical-null calibration or
  independent replication was performed. Small candidate sets limit power.
- GO significance neither controls the nominal pair-level false discoveries nor
  validates individual pairs, ligand–receptor signaling, or causal mediation.
  A negative GO result is not evidence that all candidates are false positives.

## GO implementation

Ensembl IDs map to all Entrez mappings, deduplicated at the Entrez level;
unmapped candidates and ontology coverage are recorded per screen. No candidate
has multiple Entrez mappings. The combined backgrounds contain 19 (eligible)
or 223 (all-tested) Ensembl IDs with multiple mappings; those mappings are retained.

Each ontology uses its annotated background and candidates. One-sided
hypergeometric tests use inherited GO annotations and term sizes of 10–500
background genes. The primary BH family contains all query-overlapping,
size-eligible terms within one screen/ontology, matching the DOSE calculation
used by ERC. Additional columns provide BH over all size-eligible terms and BH
pooled across 15 query-overlap families. Corrections are separate for the two
background analyses. The Python adaptation uses BH only, not the additional
Storey q-value cutoff in `enrichGO`; it is not an exact software replication.

Plotting occurs after correction: BH p < .05, at least two candidate genes,
and the union of each screen's top three terms (ranked by adjusted p, raw p,
then GO ID). Other qualifying screen associations for those terms are shown.
GO hierarchy redundancy is retained and must not be interpreted as independent
pathway support. Summary cells show minimum adjusted p, not effect size.

Reference scripts: `LFF_spatial_ERC/code/22_Mediation/05_Mediation_general_plots.R`
and `06_Mediation_GO.R`, inspected at commit
`9d8d75d9df7e3052fd65ca558584ff7123cb3dfa`; DOSE 4.4.0 `enricher_internal`.
Annotations: [org.Hs.eg.db 3.22.0](https://bioconductor.org/packages/3.22/data/annotation/html/org.Hs.eg.db.html)
and [GO.db 3.22.0](https://bioconductor.org/packages/3.22/data/annotation/html/GO.db.html),
GO source date 2025-07-22. Database metadata, hashes, software versions and compact
annotation snapshots accompany both results. All downstream calculations ran locally.

## Deliverables

Under `plots/18_mediation/nominal_p05/`, PDF and SVG are vector outputs; PNG is 300 dpi:

- `figure_6B_gene_overlap`: outcome-gene UpSet plot.
- `figure_6C_background_comparison`: side-by-side GO summary for both backgrounds.
- `figure_6C_GO_summary` and `figure_6C_GO_{CC,BP,MF}`: eligible-DEG results.
- `all_tested/figure_6C_GO_summary` and `all_tested/figure_6C_GO_{CC,BP,MF}`:
  all-tested results, including thresholded BP and MF dotplots.

Data under `processed-data/18_mediation/run_2026-10-08_limma/`:
`downstream_nominal/` and `downstream_nominal_all_tested/` contain complete GO
results, significant results, mapping/universe/annotation tables, overlap tables
and provenance. `report/GO_background_comparison.tsv` summarizes both analyses.
`GO_eligible_universes.tsv` lists the selected background in each directory;
`background_mode` in provenance and group summaries identifies which one.
The committed `primary/all_pairs.tsv.gz` contains aggregate fitted statistics
needed to reproduce the selection and backgrounds; it contains no donor-level rows.
Its archived `screen_hit` flags describe the previous FDR screen, not the nominal
selection, which is recomputed from the four gates above.

## Reproduction and verification

Install `requirements-nominal-figures.txt` in a Python environment. From the
project root:

```bash
python code/analysis/18_mediation/export_nominal_results.py \
  --source-run processed-data/18_mediation/run_2026-10-08_limma
python code/analysis/18_mediation/downstream_nominal.py \
  --source-run processed-data/18_mediation/run_2026-10-08_limma \
  --annotation-snapshot processed-data/18_mediation/run_2026-10-08_limma/downstream_nominal \
  --plots plots/18_mediation/nominal_p05
python code/analysis/18_mediation/downstream_nominal.py \
  --source-run processed-data/18_mediation/run_2026-10-08_limma \
  --annotation-snapshot processed-data/18_mediation/run_2026-10-08_limma/downstream_nominal_all_tested \
  --background all_tested --plots plots/18_mediation/nominal_p05/all_tested \
  --compare-backgrounds
Rscript code/analysis/18_mediation/tests/verify_downstream.R \
  processed-data/18_mediation/run_2026-10-08_limma eligible
Rscript code/analysis/18_mediation/tests/verify_downstream.R \
  processed-data/18_mediation/run_2026-10-08_limma all_tested
```

Initial annotation extraction uses `--orgdb PATH/org.Hs.eg.sqlite` and
`--godb PATH/GO.sqlite` instead of `--annotation-snapshot`; databases are opened
read-only. Saved snapshots avoid installing Bioconductor for reproduction.

Independent base-R checks passed for all 7,703 reported term tests: input gates,
background membership, mapped gene sets, ontology denominators, term membership,
hypergeometric probabilities, within-family/expanded/pooled BH corrections and
summary counts. The original eligible-background probabilities reproduce exactly.
The existing fitted-output verifier also passed for five screens and 44,247
pairs, including full-family FDR, sample counts and shared baselines. These checks
verify implementation and saved results; they do not establish model assumptions
or causal validity. Figures were visually checked for labels, legends and clipping.
