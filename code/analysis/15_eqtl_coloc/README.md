
# eQTL and Colocalization Analysis: DLPFC Spatial Transcriptomics

This directory contains the analysis pipeline for identifying expression quantitative trait loci (eQTLs) and performing colocalization (coloc) analysis on spatial transcriptomics data from the human DLPFC. The analysis focuses on Schizophrenia (SCZ) versus Control (NTC) donors across specific spatial domains (SpDs) and proteogenomic microenvironments (SPGs).

## 1. Data Inputs

The analysis utilizes **5 pseudobulked objects** (SpatialExperiment/SingleCellExperiment) located on JHPCE.
The path for the data is `/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/rds/`.

### A. Spatial Domain (SpD) data

In this study we have the following spatial domains: SpD01_WMtz, SpD02_L3-4, SpD03_L6, SpD04_WM, SpD05_L5, SpD06_L2-3, SpD07_L1-M where WMtz=White_Matter_transition_zone, WM=White_Matter, L=Layer, M=Meninges

* **File:** `07_dx_pseudobulk/sce_pseudo_PRECAST07_donor_spd.rds`
* **Description:** Layer-level/Domain-level pseudobulk data.
* **Action:** This object must be subsetted by Spatial Domain (SpD) before analysis. There are 63 donors.

### B. Proteogenomics microenvironment (SPG) data

These files represent PNN-informed domain labels and other microenvironments:

* **Neuropil:** `PB_dx_spg/pseudo_neuropil_pos_donor_spd.rds`
* **Neuronal:** `PB_dx_spg/pseudo_neun_pos_donor_spd.rds`
* **PNN:** `PB_dx_spg/pseudo_pnn_pos_donor_spd.rds`
* **Vasculature:** `PB_dx_spg/pseudo_vasc_pos_donor_spd.rds`

**Note:** Ideally, these files should contain one column per donor (N=63). If columns > 63, the data is likely split by domain and requires re-assessment regarding aggregation or running per-domain-microenvironment subsets.

## 2. Analysis Workflow

The pipeline utilizes `tensorQTL` (Python) for eQTL discovery and `coloc` (R) for colocalization with GWAS data.

### Data Preparation & Covariates

For each expression subset (e.g., *SpD01 only* or *Neuropil only*), we prepare the phenotype (expression) and covariate tables.

**Covariates to include:**

1. **Demographics:** Sex, Age.
2. **Technical:** Slide ID (for batch effect correction).
3. **Genotype:** Top 5 Ancestry Principal Components (snpPCs) derived from DNA genotype data.
4. **Expression Heterogeneity:** Gene expression PCs computed individually on each subset (number determined via `sva::num.sv()`).

### tensorQTL Models

We run two distinct eQTL models using `tensorQTL`:

**Model A: Base Model (Standard eQTL)**

* **Design:** `Expression ~ Genotype + Dx + Age + Sex + Slide + snpPCs + exprPCs`
* **Goal:** Identify cis-eQTLs adjusting for Diagnosis (SCZ/NTC).
* **Modes:** Nominal, Cis, Independent.

**Model B: Interaction Model (G x Dx)**

* **Design:** Tests for the interaction between Genotype and Diagnosis.
* **Goal:** Identify eQTLs that differ in effect size between SCZ and Control donors.
* *Note:* We anticipate limited power here, but it serves as an exploratory step.

> **Constraint:** We are *not* using Linear Mixed Models (LMM) to account for repeated measures (e.g., multiple layers per donor) because `tensorQTL` requires unique genotypes per sample. Therefore, analyses are stratified by domain(SpD) or microenvironment (SPG).

### Annotation & Colocalization

1. **GWAS Integration:** Annotate whether significant eQTL SNPs are risk SNPs in the PGC3 SCZ GWAS (European Ancestry only).
2. **Colocalization:** Run `coloc` on the results to assess overlap between eQTL signals and SCZ GWAS signals.

### Integration with DEGs

Join eQTL/coloc results with Differential Expression (DEG) results to identify "trifectas":

* Gene is a **DEG**.
* Gene is an **eQTL**.
* The eQTL SNP is a **GWAS risk SNP**.

## 3. Script Structure

* `01_prep_inputs.R`: Loads RDS files, subsets by SpD or SPG, calculates expression PCs, and exports formatted BED/Covariate files for tensorQTL.
* `eqtl_mapping.py`: python script running tensorQTL models and modes, used by all array jobs
* `02_run_tensorqtl_spd.sh`: Array job script for running tensorQTL on SpD aggregated data
* `03_run_tensorqtl_spg.sh`: Array job script for running tensorQTL on SPG aggregated data
* `04_annotate_gwas_DE.R`: Checks overlap with PGC3 SCZ GWAS variants and DE results
* `05_coloc.R`: Runs colocalization analysis on significant results.

## 4. Requirements

* **Software:** `tensorQTL` (Python), R (`SpatialExperiment`, `sva`, `coloc`).
* **Genotype Data:** donor genotype PLINK files (European Ancestry subset).
* **Donor Info:** 63 donors (mostly EA), SCZ and NTC diagnosis.
