# Ling et al. (2024) dlPFC snRNA-seq: count matrix download and SCE construction

`download_ling2024_counts.R` downloads the **per-library, donor-assigned count matrices** (87 files, ~41 GB) from:

> Ling E, Nemesh J, Goldman M, *et al.* A concerted neuron–astrocyte program declines in ageing and schizophrenia. *Nature* 627, 604–611 (2024). https://doi.org/10.1038/s41586-024-07109-5

## Running the script

Requires R (base R only, no packages) and internet access to `data.nemoarchive.org`.

```bash
Rscript code/analysis/16_integration_ling/Download_ling_et_al/00_download_ling2024_counts.R             # saves to ./ling2024_counts
Rscript code/analysis/16_integration_ling/Download_ling_et_al/00_download_ling2024_counts.R /Users/bguo6/GitHub/spatialDLPFC_SCZ/raw-data/ling_et_al   # or a folder of your choice

Rscript 00_download_ling2024_counts.R /Users/bguo6/GitHub/spatialDLPFC_SCZ/raw-data/ling_et_al   # or a folder of your choice
```

- Each file is saved to `<name>.part` first, checked against the server's file size, then renamed.
- A failed file prints `FAILED`, and the script continues with the next file.
- Re-running is safe: completed files are skipped and only missing or failed ones are downloaded again. If the network drops entirely, the script stops with an error; re-run it to resume.
- R's default 60-second download timeout is raised to 2 hours, because the files are 0.25–0.7 GB each.

## Building a SingleCellExperiment (`01_create_sce.R`)

`01_create_sce.R` builds a `SingleCellExperiment` for a **59-donor subset**: 3 batches from 2021, 30 SCZ and 29 NTC donors. It reads the level-2 h5 files of those batches (see [File levels](#file-levels-in-counts)), writes the counts to HDF5 one library at a time, and attaches per-nucleus and donor metadata.

```bash
Rscript code/analysis/16_integration_ling/Download_ling_et_al/01_create_sce.R
```

Requires `here`, `rhdf5`, `Matrix`, `HDF5Array`, `SingleCellExperiment`, `tidyverse` and `sessioninfo`. On a laptop it takes about 9 minutes and peaks at about 8.6 GB of RAM.

### Inputs (in `raw-data/ling_et_al/`)

| File | What it provides | Source |
|---|---|---|
| `<lib>.intronic_…gmg.donors.h5` | Gene × nucleus counts for one library | `00_download_ling2024_counts.R` |
| `BA46.jointMetadata.txt.gz` (~70 MB, 1.23 M rows) | One row per nucleus: `DONOR`, cell-type labels (`predClass`, `cellclass`, `subclass`, `subtype`), `doublet`, and QC metrics such as `NUM_TRANSCRIPTS`, `pct_intronic` and `frac_contamination` | `processed/other/` on NeMO. Downloaded by `01_create_sce.R` if missing. |
| `SZvillage_donorMetadata.txt` | One row per donor (191): `Study_ID`, `Donor`, `Sex`, `Age`, `PMI`, `Schizophrenia`, and 0/1 indicators `Batch 1`–`Batch 11` | `processed/other/` on NeMO |

The h5 files contain only barcodes, without donor labels, so `BA46.jointMetadata.txt.gz` is required to link each nucleus to its donor.

### Choice of batches

`batches_keep` is set to the three 2021 batches `05-21-2021`, `05-28-2021` and `09-02-2021`: **59 donors (30 SCZ, 29 NTC)**, processed within 104 days. Each batch is a pool of ~20 donors, and every library (reaction) in a batch contains all of that batch's donors, so a donor subset has to be chosen batch by batch.

Donors per batch, from the joint metadata:

| Batch | SCZ | NTC | Donors shared with other batches |
|---|---|---|---|
| 05-07-2019, 07-10-2019, 08-09-2019, 08-27-2019, 09-18-2019 | 10 | 10 | 1–5 with each of the others in this group |
| 10-16-2019 | 10 | 10 | 1 with 07-10-2019, 1 with 05-28-2021 |
| **05-21-2021** | 10 | 10 | none |
| **05-28-2021** | 10 | 9 | 1 each with 05-07-2019, 08-09-2019, 10-16-2019 |
| **09-02-2021** | 10 | 10 | none |
| 09-03-2021 | 9 | 10 | none |
| 09-22-2021 | 7 | 7 | none |

The kept batches share no donors with each other, so each donor appears in exactly one batch. The script stops if two kept batches share a donor or a batch lacks either diagnosis.

Alternatives that were considered:

| Batches | Donors (SCZ/NTC) | Days apart | Trade-off |
|---|---|---|---|
| **05-21-2021 + 05-28-2021 + 09-02-2021 (used)** | 59 (30/29) | 104 | Same year. 05-28-2021 has 9 NTC. 05-21-2021 used the stricter 0.300 nuclei cutoff. |
| 09-18-2019 + 10-16-2019 + 05-21-2021 | 60 (30/30) | 611 | The closest set with exactly 10 + 10 per batch and no shared donors; mixes 2019 and 2021 |
| 08-09-2019 + 09-18-2019 + 10-16-2019 | 59 (30/29) | 68 | Same year, same cutoff; one NTC donor (S07629) is in both 08-09 and 09-18 |

No three batches give exactly 60 different donors from a single year, because the 2019 batches re-use donors and only two 2021 batches have 10 SCZ + 10 NTC.

### What the script does

1. **Subset the metadata** to the kept batches and check that no donor is in two kept batches.
2. **Read gene and barcode names** from the 24 h5 files. Each file lists only the genes detected in that library (31,902 genes shared by all 24, 40,526 in the union). A gene missing from a library means zero counts there. Rows are gene symbols only, because the files contain no Ensembl IDs.
3. **Drop nuclei without a joint-metadata row** (10,423 of 403,516, 2.6%). They have no donor or cell-type label, so they cannot be assigned to a sample.
4. **Write the counts to HDF5.** The 24 libraries hold about 1.6 billion non-zero counts, close to the `dgCMatrix` limit of 2³¹ − 1 and too large to `cbind()` in memory. Each library is read with `rhdf5` (the files have no `features/id` or `genome`, so `DropletUtils::read10xCounts()` cannot read them), re-indexed onto the gene union, and appended to a TENx-format HDF5 file.
5. **Build `colData`**: the joint-metadata columns, then the donor-metadata columns (`Batch 1` becomes `Batch_1`), plus `dx` (`NTC`/`SCZ` factor) and `batch` (factor).

### Output

`processed-data/rds/16_integration_ling/sce_ling2024_3batches/`, with two files that must stay in the same folder:

| File | Contents |
|---|---|
| `counts.h5` (2.2 GB) | Genes × nuclei integer counts, TENx sparse HDF5 (group `counts`) |
| `sce.rds` | The SCE without its assay. `saveRDS()` on an HDF5-backed assay stores an absolute path, and a `TENxMatrix` path cannot be reset, so the counts are attached when loading. |

Load it with:

```r
sce <- readRDS(file.path(dir, "sce.rds"))
assay(sce, "counts") <- HDF5Array::TENxMatrix(file.path(dir, "counts.h5"), group = "counts")
```

| Slot | Contents |
|---|---|
| `assay(sce, "counts")` | 40,526 genes × 393,093 nuclei (`TENxMatrix`) |
| `rowData(sce)` | `gene_name` |
| `colData(sce)` | `library`, `barcode`, `batch`, `DONOR`, `dx`, `cellclass`, `subclass`, `subtype`, the other joint-metadata columns (QC metrics, `predClass`, `doublet`), then `Study_ID`, `Sex`, `Age`, `PMI`, `Schizophrenia`, `Batch_1`–`Batch_11` |
| `metadata(sce)$source` | Citation, DOI, kept batches and the input file names |

### Notes and limitations

- **Nuclei without metadata.** The server-side README says the level-2 files hold nuclei selected after CellBender and filtered to barcodes confidently assigned to a single donor, yet about 2.5% of nuclei have no joint-metadata row. Most likely they were removed in later QC, but this is unconfirmed. They are dropped here.
- **Outlier donors are not removed.** See [Caveats for analysis](#caveats-for-analysis).
- **Nuclei without a cell class.** About 0.9% of all joint-metadata rows (all batches) have `cellclass = NA`. They are kept in the SCE and dropped by `02_celltype_dx_DE_limma-voom.R`.

## Cell-type DE (`02_celltype_dx_DE_limma-voom.R`)

SCZ vs NTC differential expression per cell type in the 60-donor SCE, with pseudobulk limma-voom.

```bash
Rscript code/analysis/16_integration_ling/Download_ling_et_al/02_celltype_dx_DE_limma-voom.R
```

Requires `here`, `HDF5Array`, `SingleCellExperiment`, `edgeR`, `limma`, `tidyverse` and `sessioninfo`.

- **Pseudobulk** counts are summed per donor × library (reaction) × cell type. Each donor is pooled into all 8 reactions of its batch, so one donor gives up to 8 pseudobulk samples per cell type. Samples with fewer than 10 nuclei are dropped. The first aggregation, at subclass level, reads all of `counts.h5` and is cached as `pb_donor_library_subclass.rds`.
- **Cell types** are tested at two levels: `cellclass` (7 major classes) and `subclass` (the neuronal and astrocyte subclasses; the other classes have none). A cell type is skipped if it has fewer than 5 donors in either group.
- **Model:** `~ dx + Age + Sex + batch`. PMI is missing for 5 of the 191 donors and is left out. Genes are filtered with `filterByExpr(group = dx)` and normalized with TMM.
- **Donor random effect:** `duplicateCorrelation(block = DONOR)` with two rounds of voom, as in the limma user guide. The SCZ effect is tested at the donor level, not the reaction level.
- **Significance:** BH FDR < 0.05 on the `dxSCZ` coefficient.

Outputs:

| File | Contents |
|---|---|
| `processed-data/rds/16_integration_ling/dx_DE_3batches/limma_voom_fits.rds` | Per level and cell type: `fit` (`MArrayLM`), `top` (full `topTable`), `info` |
| `…/dx_DE_all_genes.csv.gz` | All tested genes for all cell types |
| `…/dx_DE_n_sig_summary.csv` | Per cell type: DEGs up/down at FDR < 0.05, genes with nominal p < 0.05, genes tested, donors, pseudobulk samples, consensus donor correlation |
| `plots/16_integration_ling/dx_DE_3batches/n_DEG_per_celltype.{pdf,png}` | Number of DEGs up (red) and down (blue) in SCZ per cell type |

## Where the data live

| | |
|---|---|
| Directory (HTTP) | https://data.nemoarchive.org/other/grant/broad/mccarroll/transcriptome/sncell/10x_v3.1/human/processed/counts/ |
| NeMO collection | `nemo:dat-djf3v1i` ("Analysis_Open"), open access: https://assets.nemoarchive.org/collection/nemo:dat-djf3v1i |
| Parent record (cited in the paper) | `nemo:dat-bmx7s1t`: https://assets.nemoarchive.org/dat-bmx7s1t |
| Server-side README | `counts/README_modified_2024-11-26_counts.txt` (5.7 KB). Read this to confirm the details below. |

`counts/` holds 2,462 of the 4,577 files (~220 of 557 GB) in `processed/`; the other files sit in sibling folders. The raw FASTQs and the alignments (BAMs) are separate controlled-access collections (`nemo:dat-kqw2ekx`, `nemo:dat-e0svrgv`).

## Libraries and batches

The data are 87 10x Chromium 3′ v3.1 libraries ("reactions") from dlPFC (BA46). They come from 11 batches, and each batch is a pooled "village" of ~20 donors that were demultiplexed by SNP genotype. Library names follow `BA46_<batch date>_rxn<N>`.

| Batch | Reactions | Intronic cutoff in file name | Note |
|---|---|---|---|
| 05-07-2019 | rxn2–8 (7) | 0.200 | names end in `_MyelinRemoved` |
| 07-10-2019, 08-09-2019, 08-27-2019, 09-18-2019 | rxn1–8 | 0.200 | |
| 10-16-2019 | rxn1–8 | 0.200 | batch used to train the paper's cell-type classifier |
| 05-28-2021, 09-02-2021, 09-22-2021 | rxn1–8 | 0.200 | |
| 05-21-2021, 09-03-2021 | rxn1–8 | **0.300** | stricter nuclei cutoff |

The directory also contains `2019-03-01_BA32_NucleiDropulationPilot_*`, a pilot from a different region (BA32) that is not part of the paper. The script does not download it.

## File levels in `counts/`

The descriptions below are inferred from the file names and the paper's Methods. The server-side README has not been checked yet.

| Level | File pattern | Files | Size | What it contains |
|---|---|---|---|---|
| 1. Full library | `BA46_<lib>.h5` | 87 | ~100 GB | All barcodes for a library before nuclei are selected, so it likely includes empty droplets. Least processed. |
| **2. Donor-assigned nuclei** (downloaded by this script) | `<lib>.intronic_0.<cut>-1.000.10X.exonic+intronic.gmg.donors.h5` | 87 | ~41 GB | Barcodes kept as nuclei, apparently those with 20% (or 30%) to 100% intronic reads. Counts are exonic + intronic, and each nucleus is assigned to a donor. Probably still contains cross-donor doublets and all donors, including later-excluded ones. |
| 3. Per cell type | `<lib>…donors.digital_expression.<celltype>.<subtype>.counts.h5` | 2,257 (26 per library) | ~77 GB | Level 2 split by the 7 major cell types and their neuronal/astrocyte subtypes. Most likely the matrices after heterotypic-doublet clusters were removed. Five rare-subtype files are absent (MEIS2 ×4, L5/6 NP ×1), presumably because those libraries had no such nuclei. |
| 4. Donor pseudobulk | `BA46.<celltype>.<subtype>.noOutliers.counts.metacells.h5` | 26 | ~1.2 GB | Donor × gene summed counts across all batches, with the paper's 11 outlier donors removed. Closest to the input of the paper's latent-factor analysis. |

The 26 cell-type groups are:

- **Astrocyte:** All, Fibrous, Interlaminar, Protoplasmic
- **Endothelia:** All
- **GABAergic:** All, LAMP5, MEIS2, PVALB, SNCG, SST, SSTCHODL, VIP
- **Glutamatergic:** All, L23IT, L4IT, L56NP, L5ET, L5IT, L6CT, L6IT, L6ITCar3, L6b
- **Microglia:** All
- **Oligodendrocyte:** All
- **Polydendrocyte:** All

## Caveats for analysis

- **Outlier donors.** The paper excluded 11 of 191 donors (4 control, 7 schizophrenia). Only the level-4 files are labelled `noOutliers`, so remove these donors from level-2 and level-3 data yourself (Supplementary Table 1 and Supplementary Fig. 1d of the paper).
- **Replicate donors.** Some donors were deliberately processed in more than one batch for quality control, so the same donor ID can appear in several libraries.
- **Points to confirm in the server-side README:** whether CellBender ambient-RNA removal is applied at level 1 or 2, whether level 3 is doublet-filtered, the meaning of `gmg`, and the internal HDF5 layout of the files.