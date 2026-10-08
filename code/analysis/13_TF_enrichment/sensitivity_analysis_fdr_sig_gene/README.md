Sensitivity analysis of the ChEA3 TF enrichment in `../`, using layer-restricted DEGs with FDR < 0.10 (`adj.P.Val < 0.10`) instead of nominal p-value < 0.05.

Steps:
1. run 00_prepare_files_for_ChEA3.R to create up-/down-regulated DEGs (FDR < 0.10) per spatial domain.
  - Domain-direction gene sets with fewer than 5 DEGs are skipped (`min_n_genes`), so only a subset of domains are included.
2. run 01_run_ChEA3_API.R to run ChEA3 (Version 3) through its API, instead of the manual portal step in the main analysis.
3. run 10_prepare_ChEA3_output.R and 11_viz_up_TF.R to compile and visualize the TF results.

Outputs:
- `processed-data/rds/13_TF_enrichment/sensitivity_analysis_fdr_sig_gene/`
- `plots/13_TF_enrichment/sensitivity_analysis_fdr_sig_gene/`
