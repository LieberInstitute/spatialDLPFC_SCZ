## Sparsity-aware data preparation strategy for pseudobulk eQTL mapping with tensorQTL

Pseudobulked data are already available in the provided RDS files. However for eQTL mapping the expression matrix (`counts`) must be filtered, normalized and transformed appropriately using different approaches for the two ways it is used in tensorQTL:

  *  expression BED file: TMM normalization (with log transform CPM) followed by inverse normal transformed (INT)
  *  covariate file (expression PCs): TMM normalization (with log transform CPM) then Z-score standardization

Xue et al. paper (Genome Biol. 2023 Feb 23;24:33. doi: 10.1186/s13059-023-02873-5) argues for the Z-score transformation (SD=1, mean 0) for computing expression PCs, and performing PCA only on the high-variance genes (HVGs) to get the expression PCs This approach seems to also minimize the correlation between expression PCs and known covariates.
