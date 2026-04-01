# coloc-gene-overlaps

This note documents the external schizophrenia-related comparison sets used in
`07_coloc_final_explore.Rmd` for the current **SCZD_PNN** DLPFC spatial
transcriptomics colocalization analysis.

## Main framing

These comparison sets are not equivalent in tissue, developmental stage, donor
composition, or evidence type.

Use them as:

- external support
- cross-study concordance
- context for novelty

Do not describe refs 85-87 as direct replication datasets.

## Primary comparison logic

The notebook now uses two layers of comparison:

1. an **adult-focused primary comparison**
   - `SCZD_PNN`
   - `habenulaPilot`
   - `Zeng2022_ref85_SCZ`
   - `Wang2018_ref87_INT18`

2. a **broad external-context comparison**
   - `SCZD_PNN`
   - `habenulaPilot`
   - `Zeng2022_ref85_SCZ`
   - `Wen2024_ref86`
   - `Wang2018_ref87_INT18`

`Wen2024_ref86` stays in the broader context because it is a fetal cortical
reference rather than an adult-brain matched comparator.

## Ref85 update

### What changed

The old file `reference85_scz_bd_20genes_filtered.csv` is no longer the primary
Ref85 input for overlap analysis.

The notebook now uses:

- `reference85_scz_only_derived.csv`

### Derivation rule

This derived file was built from the public BREMA results table associated with
Zeng et al. 2022 using the **SCZ-present** rule:

- trait must equal `Schizophrenia`
- `Shared causal prob > 0.01`
- genes are retained even if they also appear for `Bipolar disorder` or
  `Schizophrenia/Bipolar disorder`

This means the resulting Ref85 set is best described as:

- **adult multi-region brain SCZ-focused colocalization candidates from Zeng et al.**

and not as:

- legacy mixed psychiatric panel for the primary notebook comparison
- SCZ-exclusive subset

### Derived Ref85 genes

`reference85_scz_only_derived.csv` contains 15 genes:

- `ZNF823`
- `THOC7`
- `FURIN`
- `FAM134A`
- `ZFAND2B`
- `CACNA1C`
- `CNPPD1`
- `ABCB6`
- `INO80E`
- `PCCB`
- `CNTN4`
- `RERE`
- `CLCN3`
- `B3GAT1`
- `GATAD2A`

### Provenance

Ref85 derivation source:

- paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC8852232/
- PDF: https://gwern.net/doc/genetics/heritable/2022-zeng.pdf
- BREMA results table: https://hoffmg01.dmz.hpc.mssm.edu/brema/v2/brain_related_traits.html

## External study summaries

### Current study: `SCZD_PNN`

- adult postmortem human DLPFC spatial transcriptomics
- donor-level pseudobulk across 63 donors
- current manuscript colocalization set

### `habenulaPilot`

- adult postmortem human habenula bulk RNA-seq colocalization set
- 68 samples in the local comparison table
- closest prior-study comparator in this notebook

### Ref85: `Zeng2022_ref85_SCZ`

- Zeng et al. 2022
- adult multi-region brain eQTL/GWAS colocalization candidates
- source cohorts span PsychENCODE, ROSMAP, and GTEx brain data
- 3,983 RNA-seq samples from 2,119 donors
- multi-ancestry, including 474 non-European individuals
- use as adult-brain external support, not DLPFC-matched replication

### Ref86: `Wen2024_ref86`

- developing human brain xQTL/colocalization resource
- fetal cortical RNA-seq across five cohorts
- 672 RNA-seq samples
- cross-ancestry; approximately 45% African ancestry and 42% European ancestry
- use as developmental external support, not adult-brain replication

Source:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC10029021/

### Ref87: `Wang2018_ref87_INT18`

- Wang et al. 2018 / PsychENCODE
- broad integrative SCZ risk-gene prioritization framework
- not a coloc-only list
- resource spans adult PFC, temporal cortex, and cerebellum
- use as broad SCZ-risk context

Source:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC6413328/

### Ref87 sensitivity: `Wang2018_ref87_INT17`

- high-confidence subset of `INT18`
- sensitivity only
- should not replace `INT18` as the broad primary Ref87 comparator

## Interpretation guidance

Preferred language:

- external support
- cross-study concordance
- adult-focused context
- developmental context
- integrative SCZ-risk context

Avoid language implying:

- direct replication
- same-evidence benchmarking
- one-to-one tissue matching across all sets
