# Cross-study Gene-List Overlap Inputs

This note documents the datasets currently used by
`07_coloc_final_explore.Rmd` for the cross-study gene-symbol overlap analysis.
It is intended as a factual record of the comparison inputs and how they are
used in the current code.

## Current implementation in code

The current overlap code uses:

- a primary 5-way comparison across `SCZD_PNN`, `habenulaPilot`,
  `Zeng2022_ref85`, `Wen2024_ref86`, and `Wang2018_ref87_INT18`
- a secondary 4-way comparison that removes `habenulaPilot`

Matching is done on **gene symbols only**.

In code, all sets are normalized with `std_symbol_key()` and converted to a
deduplicated symbol set with `to_gene_set()`. This means:

- symbols are uppercased and trimmed
- blank / missing entries are removed
- repeated rows are collapsed to unique gene symbols
- variant IDs and Ensembl IDs are **not** used for the overlap itself

## Summary of the current inputs

| Set name | File used in code | Column used for overlap | Unique symbols used now | What the list represents |
| --- | --- | --- | ---: | --- |
| `SCZD_PNN` | derived in-memory from current notebook objects | `strong_dlpfc$gene_name` | 17 | current-study strong coloc genes (`PP4 > 0.8`) |
| `habenulaPilot` | `habenulaPilot_coloc.tsv` | `gene_name` | 16 | prior habenula coloc gene list, collapsed to unique gene symbols |
| `Zeng2022_ref85` | `reference85_scz_bd_20genes_filtered.csv` | `gene_name` | 20 | adult brain eQTL/GWAS candidate causal genes for SZ and/or BD |
| `Wen2024_ref86` | `reference86_fetal_brain_colocalized_genes.csv` | `gene_name` | 5 | developing-brain colocalized gene list |
| `Wang2018_ref87_INT18` | `Ref87_INT-18_SCZ_Risk_Gene_List.csv` | `sczgenenames` | 908 | broad PsychENCODE SCZ risk-gene prioritization list |

## 1. Current study: `SCZD_PNN`

- Study type: adult human dlPFC spatial transcriptomics with pseudobulked
  donor-level eQTL and colocalization analysis
- Brain region: dorsolateral prefrontal cortex
- Spatial contexts represented in the coloc analysis: `spd01` to `spd07`,
  `neun`, `neuropil`, `pnn`, `vasc`
- Donors currently used in this repository: 63
- Diagnosis composition: 32 `SCZ`, 31 `NTC`
- Sex: 34 male, 29 female
- Age (years): median 48.58, mean 47.03, range 23.16 to 61.39
- Gene set used for overlap: the deduplicated symbol set from current
  strong-coloc rows in `strong_dlpfc`, defined by `PP4 > 0.8`
- Evidence type of the overlap genes: **current-study strong colocalization
  genes**, not generic eGenes

Current symbols used in the overlap:

- `AC004148.1`
- `ACYP2`
- `APC2`
- `ARL17B`
- `CDHR1`
- `CTDSPL2`
- `DNPH1`
- `EFEMP1`
- `ELOVL1`
- `GOLGA6L9`
- `KANSL1-AS1`
- `LY6H`
- `RNASEH2C`
- `RPS17`
- `SETDB2`
- `SLC25A27`
- `SMG6`

## 2. Habenula comparison set: `habenulaPilot`

- File used in code: `habenulaPilot_coloc.tsv`
- Current file structure: 260 gene-variant rows, collapsed in code to 16 unique
  gene symbols for overlap
- Source data in the local habenula project: `/processed-data/18_coloc/supp_table.csv`
- Brain region: habenula
- Donors in the local habenula colData used for project context: 68
- Diagnosis composition: 35 `Schizo`, 33 `Control`
- Sex: all male
- Race: all `CAUC`
- Age at death (years): median 43.85, mean 44.14, range 20.22 to 68.00
- Relationship to current study: this is the only comparison set here with
  known donor overlap based on `BrNum` IDs; the current notebook identifies 10
  overlapping donors
- Gene set used for overlap: unique normalized symbols from the `gene_name`
  column of `habenulaPilot_coloc.tsv`
- Evidence type of the overlap genes: **habenula coloc genes**, not generic
  eGenes

Current symbols used in the overlap:

- `ANKRD45`
- `ATXN7`
- `C5ORF63`
- `CCDC122`
- `CD40`
- `DNAH10OS`
- `GABBR2`
- `IFT52`
- `LRRC37A2`
- `PABPC1L`
- `PCCB`
- `RGS16`
- `RP11-166B2.1`
- `SLC25A27`
- `TDRD6`
- `UPF1`

## 3. Ref 85: `Zeng2022_ref85`

- Reference: Zeng et al., 2022
- DOI: `10.1038/s41588-021-00987-9`
- Title: *Multi-ancestry eQTL meta-analysis of human brain identifies candidate
  causal variants for brain-related traits*
- File used in code: `reference85_scz_bd_20genes_filtered.csv`
- Column used for overlap: `gene_name`
- Current overlap set size: 20 unique normalized symbols

Study design and source material:

- adult postmortem multi-region brain eQTL meta-analysis
- 3,983 RNA-seq samples from 2,119 donors
- includes 474 non-European individuals
- combines dorsolateral prefrontal cortex samples from PsychENCODE and ROSMAP
  with 13 GTEx brain regions
- GTEx brain tissues relevant to this paper are the standard v8 brain panel:
  amygdala, anterior cingulate cortex (BA24), caudate, cerebellar hemisphere,
  cerebellum, cortex, frontal cortex (BA9), hippocampus, hypothalamus, nucleus
  accumbens, putamen, spinal cord (cervical C1), and substantia nigra

What the current comparison file represents:

- the code does **not** use all eGenes from Zeng et al.
- the code uses a local 20-gene file of **GWAS-prioritized / candidate causal
  genes** for schizophrenia and/or bipolar disorder
- this should therefore be described as an **adult brain eQTL/GWAS candidate
  causal gene list**, not as a generic eGene list

Current symbols used in the overlap:

- `ABCB6`
- `B3GAT1`
- `CACNA1C`
- `CLCN3`
- `CNPPD1`
- `CNTN4`
- `EEF1A2`
- `FAM134A`
- `FURIN`
- `GATAD2A`
- `INO80E`
- `KCTD13`
- `LL22NC03-86G7.1`
- `PBX4`
- `PCCB`
- `RERE`
- `THOC7`
- `XPNPEP3`
- `ZFAND2B`
- `ZNF823`

## 4. Ref 86: `Wen2024_ref86`

- Reference: Wen et al., 2024
- DOI: `10.1126/science.adh0829`
- Title: *Cross-ancestry atlas of gene, isoform, and splicing regulation in the
  developing human brain*
- File used in code: `reference86_fetal_brain_colocalized_genes.csv`
- Column used for overlap: `gene_name`
- Current overlap set size: 5 unique normalized symbols

Study design and source material:

- developing human brain xQTL resource
- 672 distinct developing brain donors across five cohorts
- 654 samples retained with matched genotype and cortical RNA-seq data
- developmental span: 4 to 39 post-conception weeks
- ancestry composition reported in the paper: 45% European, 25%
  Latino/admixed American, 22% African-American, 8% East Asian / Southeast
  Asian
- tissue context described in the paper as **developing human neocortex** /
  cortical RNA-seq, not adult dlPFC

What the current comparison file represents:

- the code does **not** use all fetal eGenes or all fetal xQTL genes from Wen
  et al.
- the code uses a local 5-gene file of **fetal-brain colocalized genes**
- this should therefore be described as a **developing-brain colocalization
  comparison set**, not as a generic eGene list

Current symbols used in the overlap:

- `ANKRD45`
- `CCDC122`
- `DNAH10OS`
- `LRRC37A2`
- `PCCB`

## 5. Ref 87: `Wang2018_ref87_INT18`

- Reference: Wang et al., 2018
- DOI: `10.1126/science.aat8464`
- Title: *Comprehensive functional genomic resource and integrative model for
  the human brain*
- File used in code: `Ref87_INT-18_SCZ_Risk_Gene_List.csv`
- Columns in file: `sczgenenames`, `ensembl_names`
- Column used for overlap: `sczgenenames`
- Current overlap set size: 908 unique normalized symbols

Study design and source material:

- adult brain functional-genomics resource generated by PsychENCODE
- 1,866 individuals in the integrated resource
- diagnosis composition reported in the article summary: 926 PsychENCODE
  controls, 113 GTEx controls, 558 schizophrenia, 217 bipolar disorder,
  44 autism spectrum disorder, and 8 affective disorder
- adult brain regions represented in the resource: prefrontal cortex (`PFC`),
  temporal cortex (`TC`), and cerebellum (`CB`)

What the current comparison file represents:

- the code does **not** use PsychENCODE adult-brain eGenes
- the code does **not** use a pure colocalization-gene list
- the code uses the broad `INT18` **schizophrenia risk-gene prioritization**
  file from the PsychENCODE resource
- this should therefore be described as a **broad integrative SCZ risk-gene
  reference set**

Because the current overlap code matches only on `sczgenenames`, Ensembl IDs in
the file are retained for provenance in the source CSV but are not used in the
current intersection logic.

## Evidence-class summary

The five overlap sets are not the same evidence class.

- `SCZD_PNN`: current-study strong coloc genes
- `habenulaPilot`: prior coloc genes
- `Zeng2022_ref85`: adult brain eQTL/GWAS candidate causal genes
- `Wen2024_ref86`: developing-brain colocalized genes from an xQTL atlas
- `Wang2018_ref87_INT18`: broad PsychENCODE SCZ risk-gene prioritization set

Accordingly, refs 85-87 are not all the same type of adult dlPFC colocalization
resource. In the current notebook, the comparison is a **gene-level overlap
against external schizophrenia-related prioritization resources**.
