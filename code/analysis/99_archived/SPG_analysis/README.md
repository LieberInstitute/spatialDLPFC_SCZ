Files:
1. 01-spg_spot_call.R: Calling SPG positive spots (PNN+, Neuropil+, Vasculature+, and Neun+), pseudobulked data are created for each channel.
2. plot_spg_prop.r: Create box plots to show the distribution of proportion of SPG+ spots across samples/ spatial domains.
3. SPG_DE_analysis.R: test DEG in schizophrenia among SPG+ spots 
4. plot_vocano_DE.r: given an outcome file (in csv) from pseudobulked analysis, create a volcano plot for the analysis.


pnn_L2345 includes files that specifically to do pnn analysis within a subset of spatial domains

`viz_spd_pb_data.R` creates a dot plot showing which spatial domains are still included in the pseudobulked and and how many cells are in the spatial domain. It also creates spd-specific rds files for L2/3 (spd06), L3/4 (spd02), and L5 (spd05).


Misc Files:
1. neighbor/test_find_neighbors.r: Examine if `BayesSpace:::.find_neighbors` work as expected such that we cand define PNN neighbor analysis.

TO be archived:
1. spg_spot_prop_per_spd.R: previous code
2. neighbor/01-spg_spot_call_neighbors.R: the code that Boyi actually didn't finish devleoping...
