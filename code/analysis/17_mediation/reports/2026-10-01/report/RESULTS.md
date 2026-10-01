# SCZ PTN/FGF exploratory mediation screening

Generated: 2026-10-01 15:00:01.670003
Run signature: d7fb6893ce5cf7afaa6a2d627eb7359ca2bfb80766e0ec811c96cb589953be8b

| Screen | Donors | Genes tested | Mediator Dx p | Hits | Higher priority |
|---|---:|---:|---:|---:|---:|
| ptn_vasc_neuropil | 61 | 11756 | 0.5746 | 0 | 0 |
| ptn_vasc_neun | 58 | 13209 | 0.993 | 0 | 0 |
| fgf1_neuropil_vasc | 61 | 5103 | 0.2579 | 0 | 0 |
| fgf2_neuropil_vasc | 61 | 5103 | 0.3848 | 0 | 0 |
| fgf1_neun_vasc | 58 | 5158 | 0.3433 | 0 | 0 |

17 of 17 required sensitivity tasks completed.
Overlap audit completed. 0 screens reaggregated; 5 did not require reaggregation; 0 were untestable or failed.

The screening rule requires historical and matched-baseline diagnosis p < 0.05, matched mediator diagnosis p < 0.05, conditional mediator-outcome BH FDR < 0.05, and adjusted diagnosis p >= 0.05.

Higher priority additionally requires same-direction coefficient shrinkage and BH FDR < 0.05 across all five mediator-outcome testing families. This is not a mediation-discovery FDR guarantee.

All four historical mediator nominations fail gene-wide FDR 0.05. Nominal selection and testing reuse the same cohort. A p-value crossing 0.05 is not a test of coefficient change or proof of mediation.

These are cross-sectional mixed-tissue observations with repeated SpDs per donor and potentially shared spots between SPGs. Diagnosis is not randomized; temporal ordering and unmeasured confounding remain unresolved. Results support hypotheses for follow-up, not causal or functional validation.

See screen_summary.tsv, hit_sensitivity_comparison.tsv (when hits exist), primary/all_pairs.tsv.gz, historical/reconciliation.tsv, and the sensitivity/overlap status tables for complete evidence.
