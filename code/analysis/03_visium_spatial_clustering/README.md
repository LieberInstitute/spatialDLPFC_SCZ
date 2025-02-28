Order to run files:

01_PRECAST.R/sh: running PRECAST with selected marker genes from PEC study
02_post_PRECAST_int.R: extract PRECAST labels, save as an file `processed-data/rds/spatial_cluster/PRECAST/test_clus_label_df_semi_inform_k_2-16.rds`, and merge to SPE object

plot_*: make plots for panels in Figure 1 or Supplementary figures.


# NOTE
- If needing revision:  can replace all the load spe object from `processed-data/rds/01_build_spe/fnl_spe_kept_spot_only.rds`, such that the merging steps in the code chunks can be removed.