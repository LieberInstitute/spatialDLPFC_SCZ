## fixed scientific choices; paths resolve against --project-root.
config <- list(
  version = "1.0.0",
  spds = sprintf("spd%02d", c(2, 3, 5, 6, 7)),
  min_donors = 10L, min_per_dx = 5L, min_residual_df = 5L,
  min_spots = 10L, dx_p = 0.05, mediator_dx_p = 0.05, mediator_q = 0.05,
  seed = 172026L,
  pb = c(vasc = "processed-data/rds/PB_dx_spg/pseudo_vasc_pos_donor_spd.rds",
         neun = "processed-data/rds/PB_dx_spg/pseudo_neun_pos_donor_spd.rds",
         neuropil = "processed-data/rds/PB_dx_spg/pseudo_neuropil_pos_donor_spd.rds"),
  historical = c(vasc = "processed-data/spg_pb_de/test_SPD_pseudo_vasc_pos.csv",
                 neun = "code/analysis/dx_deg_spg_neun/neun-dx_DEG-GM.csv",
                 neuropil = "code/analysis/dx_deg_spg_neuropil/neuropil-dx_DEG-GM.csv"),
  donor_meta = "processed-data/ref/donor_meta.tsv.gz",
  spots = "processed-data/rds/01_build_spe/fnl_spe_kept_spots_only.rds",
  expected = data.frame(context = c("vasc", "neun", "neuropil"),
                        genes = c(5215L, 16768L, 11834L),
                        nominal = c(440L, 1789L, 1669L)),
  source_reference = "LFF_spatial_ERC/code/22_Mediation at 9d8d75d9df7e3052fd65ca558584ff7123cb3dfa"
)
