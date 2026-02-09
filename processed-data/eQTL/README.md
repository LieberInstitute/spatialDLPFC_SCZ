## tensorQTL input/output files

For each context (7 SpD + 4 microenvironments), input and output files use a specific
recognizable file name prefix.

### Input files

`tqtl_in` folder has the input expression and covariate files for tensorQTL run,
as used by the script `code/analysis/15_eqtl_coloc/eqtl_mapping.py`

Genotype data was made available in processed-data/genotypes and SNP PCs 
were generated with `code/analysis/15_eqtl_coloc/00_get_SNP_PCs.sh`

Expression-related data was prepared by `code/analysis/15_eqtl_coloc/01_prep_inputs.R`

### tensorQTL output files

These are in the `tqtl_out` folder, as follows, for each context:

  - `*.gene.cis_qtl_pairs.chr*.parquet` : `map_nominal` outputs per chromosome
  - `*.gene.eqtl-nominal_FDR05.tab.gz`, `*.gene.eqtl-nominal_p001.tab.gz` : 
       aggregated `map_nominal` outputs filtered by BH p-adjusted <.05 
       and nominal p<.001 respectively
  - `*.gene.map_cis.tab.gz` : `map_cis` permutation-based output 
  - `*.gene.map_independent.txt.gz` : `map_independent` output for FDR<0.05 map_cis features
  - `*.gene.DxINT.cis_qtl_top_assoc.txt.gz` : output for `map_nominal` with DX interaction term testing

The columns for these formats are documented here:
 https://github.com/broadinstitute/tensorqtl/blob/master/docs/outputs.md

Exploratory code for these data can be found in `code/analysis/15_eqtl_coloc/03_eqtl_explore.Rmd`,
and  `03a_nominal_eQTLs.Rmd`,`03b_eQTL_boxplots.Rmd`. 


### Colocalization output files

In the `coloc` folder, the `coloc.abf` result files are:

  - `coloc_*_gene.qs2` : unfiltered coloc.abf results per context
  - `coloc_*_gene.flat.tab.gz` : flattened files filtered by minimal PP3 and PP4 criteria, with annotations/classification codes 

Please see the `code/analysis/15_eqtl_coloc/05_coloc_explore.Rmd` for the thresholds and annotation codes used, and summary plots. 
