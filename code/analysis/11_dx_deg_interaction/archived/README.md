1. create_limma_fit_interaction.R - specify interaction model, run limma fit
1. hypo_test_dx_deg_across_spd - construct contrast matrix examining any SpD have non-zero  FC, save laminar-specific log2FC estimates
1. hypo_test_logFC_varying_across_spd.r - Test if the effect size is the same across spatial domains. **Spatially Divergent Gene**
1. Create_digital_null_gene.r - permute the psuedobulked gene with in each individual gene and calculate the logFC for each spd
1. upset_plot_layer_vulnerability.r - an upset plot showing layer-specific DEGs.



NOTE:
- Spatially divergent gene
  - fdr <= 0.05 is *HNRNPH3*, *ALDH1A1*
  - fdr <= 0.10 is *HNRNPH3*, *ALDH1A1*, *SOD2*, *CCNI*
