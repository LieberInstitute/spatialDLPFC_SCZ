# coloc-gene-overlaps

We call this current repository study the "SRT DLPFC" study, or SCZD_PNN

## Purpose

This note documents the external schizophrenia-related reference gene sets to compare against the current **SCZD_PNN** SRT DLPFC colocalization results, together with the planned **habenulaPilot** comparison set.

The immediate goals are:

1. build a **5-way Venn diagram** and using **gene names only** for:

   * current SRT DLPFC SCZ study (SCZD_PNN strong coloc genes)
   * `habenulaPilot_coloc.tsv`
   * Ref 85-derived list
   * Ref 86-derived list
   * Ref 87-derived list

Ref 85, 86, 87 are references used in the habenula pilot paper to compare the gene lists from there with the habenula coloc results, so we want to intersect our coloc results with all these datasets. 

Also, build a table with gene lists for each of the intersections between these datasets (when the intersection has lower than 20 genes). the table should also have a row with the list of genes that are unique to to current SCZD_PNN, see 2. below:

2. determine which colocalized genes are **unique to the current SRT DLPFC study** - inclde in the 

3. summarize overlap with prior datasets without overstating replication, especially because the external references differ in tissue, developmental stage, donor composition, and evidence type.

---

## Input datasets for overlap analysis

### 1. Current study (SRT DLPFC)

* **Name to use in plots/text:** `SCZD_PNN`
* **Description:** current SCZ case-control spatial transcriptomics study; use the set of colocalized genes from the present manuscript analysis
* **Input file:** to be defined locally
* **Matching field for Venn:** gene symbol only

### 2. Habenula comparison dataset

* **Name to use in plots/text:** `habenulaPilot`
* **Description:** prior habenulaPilot list of colocalized genes
* **Input file:** `habenulaPilot_coloc.tsv`
* **Matching field for Venn:** gene symbol only

Note: this is the only comparison expected to support manuscript language about prior bulk RNA-seq analyses with **shared donors**, assuming that is correct for the habenulaPilot comparison in the Methods. The three external references below are **independent** reference resources and should be treated as cross-study concordance, not same-donor replication.

---

## External reference datasets

### 3. Ref 85

* **Publication label:** **Zeng et al., 2022**
* **Full citation label:** *Multi-ancestry eQTL meta-analysis of human brain identifies candidate causal variants for brain-related traits*
* **Reference number in habenula paper:** **85**
* **File to use:** `reference85_scz_bd_20genes_filtered.csv`
* **Planned set name:** `Zeng2022_ref85`
* **Evidence type:** adult multi-region brain eQTL/GWAS colocalization candidates for **SZ/BD**, not strict SCZ-only
* **Important caution:** this file should be described as an **adult SZ/BD colocalization candidate list**, not as a pure SCZ-only list

Why this label:

* the Zeng study analyzed 3,983 RNA-seq samples from 2,119 donors across PsychENCODE, ROSMAP, and GTEx brain data, and identified candidate causal variants and genes for brain-related traits.
* for schizophrenia and bipolar disorder, the paper states that it identified **20 genes** predicted to confer risk for one or both disorders.
* this is therefore useful for descriptive overlap with current SRT colocalization hits, but it is not a clean SCZ-only benchmark.

### 4. Ref 86

* **Publication label:** **Wen et al., 2024**
* **Full citation label:** *Cross-ancestry atlas of gene, isoform, and splicing regulation in the developing human brain*
* **Reference number in habenula paper:** **86**
* **File to use:** `reference86_fetal_brain_colocalized_genes.csv`
* **Planned set name:** `Wen2024_ref86`
* **Evidence type:** developing-brain SCZ colocalization set
* **Genes in file:** `ANKRD45`, `CCDC122`, `DNAH10OS`, `LRRC37A2`, `PCCB`

Why this label:

* the Wen study is a **developing human brain** cross-ancestry xQTL atlas, spanning prenatal developmental stages rather than adult postmortem bulk brain.
* it is appropriate for descriptive overlap with SCZ-prioritized genes, but not as an adult bulk-brain matched comparator.

### 5. Ref 87

This reference should be tracked as **two related files**, but only one of them belongs in the 5-way Venn unless explicitly doing a second sensitivity plot.

#### 5a. Full Ref 87 set

* **Publication label:** **Wang et al., 2018**
* **Full citation label:** *Comprehensive functional genomic resource and integrative model for the human brain*
* **Reference number in habenula paper:** **87**
* **Official PsychENCODE filename:** `Ref87_INT-18_SCZ_Risk_Gene_List.csv`
* **Planned set name:** `Wang2018_ref87_INT18`
* **Evidence type:** broad PsychENCODE SCZ risk gene set
* **Size:** **1,111 genes**

Why this label:

* Wang et al. reported **1,111 putative SCZ-associated genes** linked to schizophrenia GWAS loci, of which **321** were a higher-confidence subset.
* this uploaded file is the correct **full Ref 87 gene set**.

#### 5b. High-confidence Ref 87 subset

* **Official PsychENCODE filename:** `Ref87_INT-17_SCZ_High_Confidence_Gene_List.csv`
* **Planned set name:** `Wang2018_ref87_INT17`
* **Evidence type:** high-confidence PsychENCODE SCZ subset
* **Size:** **321 genes**

Why keep this file:

* the official PsychENCODE resource distinguishes `Ref87_INT-18_SCZ_Risk_Gene_List.csv` from `Ref87_INT-17_SCZ_High_Confidence_Gene_List.csv`.
* this is useful as a stricter sensitivity set, but it should **not** replace the full Ref 87 set in the primary 5-way Venn if the goal is to represent Ref 87 broadly.

**Recommendation for the primary 5-way Venn:** use **`Ref87_INT-18_SCZ_Risk_Gene_List.csv`**, not `Ref87_INT-17`.

---

## Summary table

| Ref       | Publication shorthand | File                                            | Use in primary 5-way Venn? | Notes                                               |
| --------- | --------------------- | ----------------------------------------------- | -------------------------- | --------------------------------------------------- |
| Current   | SCZD_PNN     | local current-study file                        | Yes                        | current SRT SCZ coloc genes                         |
| Habenula  | habenulaPilot         | `habenulaPilot_coloc.tsv`                             | Yes                        | prior coloc list; shared-donor comparison candidate |
| 85        | Zeng2022_ref85        | `reference85_scz_bd_20genes_filtered.csv`       | Yes                        | adult SZ/BD coloc candidates; not SCZ-only          |
| 86        | Wen2024_ref86         | `reference86_fetal_brain_colocalized_genes.csv` | Yes                        | developing-brain SCZ coloc genes                    |
| 87        | Wang2018_ref87_INT18  | `Ref87_INT-18_SCZ_Risk_Gene_List.csv`                 | Yes                        | full 1,111-gene PsychENCODE SCZ risk list           |
| 87 subset | Wang2018_ref87_INT17  | `Ref87_INT-17_SCZ_High_Confidence_Gene_List.csv`      | No, secondary only         | 321-gene high-confidence subset                     |

---

## Rules for overlap analysis

### Primary Venn analysis

Use **gene symbols only** for the 5-way Venn:

* `SCZD_PNN`
* `habenulaPilot`
* `Zeng2022_ref85`
* `Wen2024_ref86`
* `Wang2018_ref87_INT18`

### QC/secondary checks

Even though the Venn should use gene names only, keep the following in mind:

* Ref 87 files also contain **Ensembl IDs**, and those should be preserved for QC and reproducibility.
* Ref 87 symbol-based matching may miss or mis-handle rows with missing or unstable symbols; therefore any final overlap counts reported in text should ideally be cross-checked against Ensembl IDs where possible.
* The previously obtained `Ref87_INT-17_SCZ_High_Confidence_Gene_List.csv` is a **subset** of `Ref87_INT-18_SCZ_Risk_Gene_List.csv`, so those two should not both appear in the same primary Venn.

---

## Interpretation constraints

These five sets do **not** all represent the same kind of evidence.

* `habenulaPilot` is the most appropriate comparator for language about prior RNA-seq analyses and potentially shared donors.
* Ref 85 is an **adult multi-region brain SZ/BD colocalization candidate list**.
* Ref 86 is a **developing-brain** SCZ-related xQTL colocalization resource.
* Ref 87 is a **PsychENCODE SCZ risk gene prioritization framework**, not a pure colocalization-only list.

Therefore, overlap with refs 85-87 should be described as:

* **cross-study concordance**
* **external support**
* **context for novelty**

and **not** as formal replication or direct evidence that bulk RNA-seq in the same setting missed the signal.

---

## Main analysis target

A main deliverable from this comparison should be:

### Genes unique to the current SRT DLPFC study

Define the **SRT-unique coloc genes** as genes present in:

* `SCZD_PNN`

but absent from:

* `habenulaPilot`
* `Zeng2022_ref85`
* `Wen2024_ref86`
* `Wang2018_ref87_INT18`

This set is the most relevant for the manuscript claim that spatial transcriptomics may reveal **context-dependent regulatory signals** not captured in prior bulk or external reference analyses.

Potential manuscript framing:

* the **habenulaPilot** comparison can support statements about signals not prioritized in earlier bulk analyses
* the comparison to refs 85-87 can support a broader statement that some SRT-prioritized genes are not represented in prior adult bulk, developmental, or integrative SCZ reference resources

---

## Planned outputs

1. **5-way Venn diagram**

   * based on gene names only
   * primary purpose: visualize overlap structure across current SRT, habenulaPilot, and refs 85-87

2. **Table of current SRT genes unique to SRT**

   * gene name
   * present in current SRT
   * absent from each of the four comparison sets
   * optional annotation column for manuscript prioritization
 Should also generate a table showing how these genes are distributed across our 11 contexts (7 SpDs + 4 SPGs):
      table should have: context, n_genes, gene_list 

3. **Optional secondary table**

   * overlap of current SRT with Ref 87 high-confidence subset (`Ref87_INT-17`)
   * useful as a stricter PsychENCODE sensitivity analysis

---

## Files currently tracked

* `habenulaPilot_coloc.tsv`
* `reference85_scz_bd_20genes_filtered.csv`
* `reference86_fetal_brain_colocalized_genes.csv`
* `Ref87_INT-18_SCZ_Risk_Gene_List.csv`
* `Ref87_INT-17_SCZ_High_Confidence_Gene_List.csv`

Once `habenulaPilot_coloc.tsv` and the current SRT DLPFC colocalized gene list are in hand, the overlap analysis can proceed.
