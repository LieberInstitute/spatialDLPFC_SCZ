## Main Files
1. create_pb_data.R - Create and save pseudobulk data for each PRECAST K
2. PB_analysis_dx_gene.R - running analysis and saving the DEG some where
3. plot_sig_gene_enrichment_PEC.R - heatmap to plot the significance level in spatialLIBD (PEC study)



## Other files
__Folders containing anscillary analysis__

*analysis_per_dx* - folder containing files that doing pseudobulk analysis for each dx group respectively.

*interaction_mdl* - analysis considering the interaction term between dx and spd to highlilght DEGs with laminar specificity.



__Files that used to explore data & Results__

* descriptive_dx_DEG_across_spd.r - summarize how many DEGs detected for each PRECAST K
* plot_PB_Genes.R - Different viz (spagetti plot and heatmap) to highlight the difference of each genes comparing  NTC and SCZ samples. 

__Analysis per dx group ___

This is the analysis that __calculate/find the spd marker genes respecitively for NTC and SCZ samples respectively__

1. create_pb_data_per_group.r
2. PB_analysis_SpD_per_dx.r
3. plot_sig_gene_enrichment_SPD_per_dx