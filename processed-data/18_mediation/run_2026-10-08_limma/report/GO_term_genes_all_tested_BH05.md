# Genes in significant GO associations

11 screen–term associations, 10 distinct GO terms and 10 distinct candidate genes. Each row lists the candidate outcome genes annotated to that term, rather than all genes annotated to the term in the background.

Selection: the existing **all-tested target-gene background** results with BH-adjusted p < .05 within each screen/ontology. Statistics are unchanged. The eligible-DEG background has no significant terms. BP = biological process; MF = molecular function.

[Association CSV](GO_term_genes_all_tested_BH05.csv) · [One gene per screen–term row](GO_gene_term_membership_all_tested_BH05.csv)

| Screen | GO term | GO ID | Candidate genes | Count | BH-adjusted p |
| --- | --- | --- | --- | ---: | ---: |
| PTN vascular to neuronal | protein serine/threonine kinase binding (MF) | GO:0120283 | **CCDC102B**, **NAIP** | 2 | 0.03101 |
| FGF1 neuropil to vascular | regulation of growth (BP) | GO:0040008 | **AGT**, **BCL6**, **FAM107A**, **GAP43**, **MT1G**, **MT3** | 6 | 0.01719 |
| FGF2 neuropil to vascular | intracellular zinc ion homeostasis (BP) | GO:0006882 | **MT1G**, **MT1M**, **MT3** | 3 | 0.00710 |
| FGF2 neuropil to vascular | cellular response to cadmium ion (BP) | GO:0071276 | **MT1G**, **MT1M**, **MT3** | 3 | 0.00710 |
| FGF2 neuropil to vascular | cellular response to zinc ion (BP) | GO:0071294 | **MT1G**, **MT1M**, **MT3** | 3 | 0.00710 |
| FGF2 neuropil to vascular | response to zinc ion (BP) | GO:0010043 | **MT1G**, **MT1M**, **MT3** | 3 | 0.00909 |
| FGF2 neuropil to vascular | regulation of growth (BP) | GO:0040008 | **AGT**, **BCL6**, **FAM107A**, **GAP43**, **MT1G**, **MT1M**, **MT3** | 7 | 0.00909 |
| FGF2 neuropil to vascular | cellular response to copper ion (BP) | GO:0071280 | **MT1G**, **MT1M**, **MT3** | 3 | 0.00909 |
| FGF2 neuropil to vascular | response to cadmium ion (BP) | GO:0046686 | **MT1G**, **MT1M**, **MT3** | 3 | 0.01312 |
| FGF2 neuropil to vascular | growth (BP) | GO:0040007 | **AGT**, **BCL6**, **FAM107A**, **GAP43**, **MT1G**, **MT1M**, **MT3**, **PLEC** | 8 | 0.01475 |
| FGF2 neuropil to vascular | response to copper ion (BP) | GO:0046688 | **MT1G**, **MT1M**, **MT3** | 3 | 0.01606 |

## Interpretation

Seven FGF2 metal-response terms contain the same **MT1G, MT1M and MT3** genes. The growth terms overlap those genes and additional candidates. Repeated annotations do not represent independent pathway evidence. These results neither establish metal exposure nor validate causal mediation.

## FGF1/FGF2 genes in the UpSet plot

The existing [exact-intersection table](../downstream_nominal/gene_overlap_intersections.tsv) lists gene symbols for each plotted column. The [membership matrix](../downstream_nominal/gene_overlap_membership.tsv) provides one gene per row and one membership indicator per screen.

The neuropil FGF1 and FGF2 screens share 16 genes: 14 in those two screens only, plus 2 also in the neuronal FGF1 screen (**FAM107A, SLC24A2**). The UpSet bars show exact intersections, so these categories are separate.

## Sources and reproduction

Source: [original significant GO table](../downstream_nominal_all_tested/GO_significant_BH05.tsv). [Methods](../../../../code/analysis/18_mediation/DOWNSTREAM.md). CSV p-values retain the original precision. Other BH columns are separate sensitivity corrections.

From the repository root:

```bash
python3 code/analysis/18_mediation/export_go_gene_tables.py \
  --source-run processed-data/18_mediation/run_2026-10-08_limma
```
