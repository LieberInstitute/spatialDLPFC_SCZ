library(SpatialExperiment)
library(spatialLIBD)

path_spe_after_spot_qc <- here::here("processed-data", "rds",
                                     "spe", "spe_after_spot_qc.rds")

# In-tissue spots only
spe <- readRDS(
  path_spe_after_spot_qc
)

vis_gene(
  spe = spe, sampleid = "V12F14-053_A1",
  geneid = "sum_umi"
  )

vis_gene(
  spe = spe, sampleid = "V12F14-053_A1",
  geneid = "sum_gene"
)


vis_grid_gene(spe,
              geneid = "sum_gene",
              pdf_file = here("plots/02_visium_qc/test_oddity",
                              "random_pattern_sum_gene.pdf") )

vis_grid_gene(spe,
              geneid = "sum_umi",
              pdf_file = here("plots/02_visium_qc/test_oddity",
                              "random_pattern_sum_umi.pdf") )

vis_grid_gene(spe,
              geneid = "expr_chrM",
              pdf_file = here("plots/02_visium_qc/test_oddity",
                              "random_pattern_expr_chrM.pdf") )

vis_grid_gene(spe,
              geneid = "expr_chrM_ratio",
              pdf_file = here("plots/02_visium_qc/test_oddity",
                              "random_pattern_expr_chrM_ratio.pdf") )

vis_gene(
  spe = spe, sampleid = "V12F14-053_A1",
  geneid = "expr_chrM_ratio"
)



