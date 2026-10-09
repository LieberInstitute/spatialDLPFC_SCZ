# PTN/FGF cross-microenvironment mediation screen: report

Run: 2026-10-08, signature `50f8ff72...`, local (R 4.5.2, limma 3.66.0, edgeR 4.8.2).
Outputs: `processed-data/18_mediation/run_2026-10-08_limma/`.

## Goal

Fig 6F nominates SCZ-linked ligand-receptor (LR) axes (PTN-PTPR, FGF-FGFR1) by gene
overlap with microenvironment-restricted SCZ-DEGs (nominal p < .05). A reviewer
asked for evidence beyond this annotation. The screen tests whether SCZ-associated
PTN/FGF expression in one SPG microenvironment statistically accounts for part of the
SCZ association of genes in another microenvironment of the same donor and SpD:

- Dx -> vascular PTN -> neuropil or neuronal genes
- Dx -> neuropil FGF1 / FGF2, or neuronal FGF1 -> vascular genes

## Design

Baron & Kenny screen adapted from the ERC/LC framework (`LFF_spatial_ERC/code/22_Mediation`),
using the manuscript SPG DE model (limma on donor x SpD logcounts, `~ Dx + age + sex + SpD`,
donor blocking) on gray-matter SpDs with >= 10 spots in both microenvironments
(272 observations / 61 donors for vascular-neuropil, 220 / 58 for vascular-neuronal).

X = diagnosis (SCZ vs NTC). M = expression of PTN, FGF1 or FGF2 in the source
microenvironment, one value per donor x SpD. Y = each gene in the target
microenvironment of the same donor x SpD, one model per gene.

| Screen | M | Y |
|---|---|---|
| PTN vascular -> neuropil | PTN, vascular | neuropil genes |
| PTN vascular -> neuronal | PTN, vascular | neuronal genes |
| FGF1 neuropil -> vascular | FGF1, neuropil | vascular genes |
| FGF2 neuropil -> vascular | FGF2, neuropil | vascular genes |
| FGF1 neuronal -> vascular | FGF1, neuronal | vascular genes |

| Step | Criterion | Threshold |
|---|---|---|
| 1. X -> Y | Y is a manuscript DEG of its microenvironment, and matched-sample c p < .05 | nominal p < .05 (FDR < .10 would leave 352 neuropil, 90 neuronal, 19 vascular targets) |
| 2. X -> M | M is a manuscript DEG in its source microenvironment (PTN vascular p = .013, FGF1 neuropil .040, FGF2 neuropil .033, FGF1 neuronal .047; FDR .31-.46) | nominal p < .05 |
| 3a. M -> Y given X | b BH FDR over all target genes in the screen | FDR < .10 |
| 3b. X -> Y given M | c' no longer significant | p >= .05 (mirrors step 1) |

The manuscript DE tables were reproduced exactly (max |dp| < 5e-13). The 2026-10-01
run (voom engine plus a re-test of step 2 on matched samples, which no mediator
passed) is superseded.

## Results

| Screen | Step 3 eligible targets | b FDR < .10 | Hits | Median c shrinkage [range] |
|---|---:|---:|---:|---|
| PTN vascular -> neuropil | 1,187 | 400 | 35 | 17% [11-37%] |
| PTN vascular -> neuronal | 1,293 | 224 | 16 | 9% [8-16%] |
| FGF1 neuropil -> vascular | 258 | 48 | 12 | 16% [12-23%] |
| FGF2 neuropil -> vascular | 258 | 119 | 26 | 19% [10-37%] |
| FGF1 neuronal -> vascular | 200 | 13 | 3 | 14% [13-14%] |

- 92 pairs meet the rule; 87 also pass b FDR < .10 pooled over the five screens.
  All show same-sign attenuation of c with a x b consistent with c.
- Attenuation is partial: median 16%, 7 pairs >= 25%. c' p values are 0.051-0.16,
  i.e. the loss of significance is marginal for most pairs.
- Threshold dependence (all rows: matched c p < .05 and c' p >= .05):

  | Step 1: Y selection | Step 3a: b | PTN -> neuropil | PTN -> neuronal | FGF1 neuropil -> vasc | FGF2 neuropil -> vasc | FGF1 neuronal -> vasc | Total |
  |---|---|---:|---:|---:|---:|---:|---:|
  | nominal p < .05 | FDR < .05 | 30 | 11 | 9 | 22 | 3 | 75 |
  | nominal p < .05 | FDR < .10 (primary) | 35 | 16 | 12 | 26 | 3 | 92 |
  | nominal p < .05 | p < .05 | 39 | 22 | 17 | 26 | 6 | 110 |
  | FDR < .10 | FDR < .05 or < .10 | 1 | 0 | 0 | 0 | 0 | 1 (HSPD1) |
  | FDR < .10 | p < .05 | 1 | 0 | 0 | 0 | 1 | 2 (+ APOLD1) |
  | FDR < .05 | any of the above | 0 | 0 | 0 | 0 | 0 | 0 |

  Step 1 drives the result: strongly DE targets keep a significant c' after
  adjustment, so the pairs come from nominally DE targets where modest attenuation
  moves p just past .05. Among eligible targets, b p < .05 occurs for 433, 368, 77,
  137 and 38 genes (screen order above) against ~59, 65, 13, 13 and 10 expected by
  chance; FDR < .10 already captures most of these associations.
- Top pairs (b FDR):
  - PTN -> neuropil: DIO2, IGSF8, KCNJ16, PHYHIP, CCND3, ARAP2, TCEAL4, CAMK1G
  - PTN -> neuronal: CAMK1G, BOD1L1, PDP1, CCDC102B, LY6H, SMARCA2
  - FGF1/FGF2 neuropil -> vascular: AGT, MT3, FAM107A, MT1M, MT1G, HINT1, SLC14A1, CHGB
  - FGF1 neuronal -> vascular: FAM107A, RGS4, SLC24A2
- The FGF -> vascular hits are dominated by astrocyte-enriched genes (AGT, MT3,
  MT1M/G, FAM107A, SLC14A1); FGF2 and PTN are themselves astrocyte-expressed in
  adult cortex. The data cannot separate PTN/FGF-mediated effects from a shared
  astrocyte abundance or state present in both microenvironments of a donor/SpD.

### Mediator diagnosis effects in the matched samples (reported, not gated)

| Mediator | Manuscript (own samples) | Matched GM samples |
|---|---|---|
| PTN vascular (neuropil set) | -0.140, p = .013 (7 SpDs) | -0.118, p = .068 |
| PTN vascular (neuronal set) | same | -0.061, p = .40 |
| FGF1 neuropil | +0.113, p = .040 | +0.089, p = .12 |
| FGF2 neuropil | +0.195, p = .033 | +0.149, p = .13 |
| FGF1 neuronal | +0.226, p = .047 | +0.143, p = .21 |

Directions match the nominations; significance is lost with the smaller matched sets.
Vascular PTN by SpD (`ptn_vasc_by_spd.tsv`, all 7 SpDs): lower in SCZ in 5 of 7
domains (L1/M -0.29, p = .015; WMtz -0.25; WM -0.22; L5 -0.21; L6 -0.16), higher in
L2/3 (+0.10) and L3/4 (+0.02); WM mean -0.23 (p = .028), GM mean -0.11 (p = .087),
WM - GM difference p = .27. The PTN signal is not WM-specific.

### Shared spots

SPG labels are not exclusive: within matched samples 28% of vascular spots are also
neuropil spots and 21% are also neuronal spots. Removing shared spots from both
pseudobulks and refitting (263 / 60 and 209 / 57 observations / donors):

- b is essentially unchanged (median refit/primary ratio 0.90-1.01; FGF1 neuronal 1.38);
  80 of 92 primary hits keep b FDR < .10. The mediator-outcome associations are not
  produced by shared spots.
- Hit membership is unstable: 33 of 92 primary hits remain hits, and the refit yields
  103 hits overall, because c' p values sit near the .05 boundary.

## Interpretation and limits

- Steps 1 and 2 both rest on nominal (p < .05) DEG selection, the threshold the
  reviewers questioned; step 2 does not replicate at p < .05 in the matched samples,
  and with FDR-selected targets only one pair remains.
- "Loss of significance" is not a test of attenuation; no indirect effect or its
  uncertainty is estimated, and attenuation is modest.
- Observational, cross-sectional data: no temporal order, and common causes
  (cell composition, astrocyte state, medication, tissue quality) can induce both
  the a and b associations. 58-61 independent donors.
- The screen supports, at a nominal level, co-variation of PTN/FGF expression with
  SCZ-associated genes across microenvironments, most plausibly through a shared
  astrocyte component; it does not demonstrate ligand-receptor signaling or causation.
