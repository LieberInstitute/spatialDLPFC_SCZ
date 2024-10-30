suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
  library(spatialLIBD)
})

# Load Data ----
## Load SPG spe object ----
spe <- readRDS(here("processed-data/rds/02_visium_qc","qc_spe_w_spg_N63.rds"))

## Load SpD data ----
finalized_spd <- readRDS(here("processed-data/rds/spatial_cluster", "PRECAST","test_clus_label_df_semi_inform_k_2-16.rds"))

## Attach SpD label to spe ----
col_data_df <- colData(spe) |>
  data.frame() |>
  left_join(
    finalized_spd,
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

# Call SPG spots ----
spe$pnn_pos <- ifelse(spe$spg_PWFA > 0.05, TRUE, FALSE)
# NOTE: neuropil spot are spots doesn't have DAPI staining
spe$neuropil_pos <- ifelse(
  spe$spg_PDAPI > 0.05 & spe$spg_PDAPI < 0.5,
  FALSE, TRUE
)
spe$neun_pos <- ifelse(
  spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3,
  TRUE, FALSE
)
spe$vasc_pos <- ifelse(
  spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,
  TRUE, FALSE
)

spe_ntc <- spe[, colData(spe)$dx == "ntc"]
head(colData(spe_ntc))


 ## Per sample & domain ----
neuropil_pseudo <-
  registration_pseudobulk(
    spe_ntc,
    var_registration = "neuropil_pos",
    var_sample_id = "sample_id",
    covars = c("age", "sex", "lot_num", "slide_id"),
    min_ncells = 10,
    pseudobulk_rds_file = here("processed-data", "image_processing", "enrichment", "neuropil_pseudo.rds")
    )

dim(neuropil_pseudo)
#[1] 15719    62	

dx_res <- registration_stats_enrichment(
  neuropil_pseudo,
  block_cor = NaN,
  covars = c("age", "sex"),
  var_registration = "neuropil_pos",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)




library(readxl)

df <- read.table(here("raw-data/images/SPG_Spot_Valid/Neuropil/maddy.csv"), header = TRUE)
 unique(df$cluster)

[1] "Synapse_ExPFC" "Synapse_In"    "LO-synapse"    "N-synapse"    
[5] "ODCjunction"   "ASCjunction"

df_sub = df[df$cluster=="LO-synapse", ]

df_neuropil <- merge(dx_res, df_sub, by = "gene")
Cneuropil = cor(df_neuropil$logFC_TRUE,df_neuropil$avg_log2FC)

png(here("plots", "image_processing", "enrichment", "Neuropil.png"), width = 800, height = 600)
plot(df_neuron$logFC_TRUE, df_neuron$logFC, 
     main = paste0("neuropil ",dim(df_neuropil)[[1]],",", Cneuropil),
     xlab = "logFC_TRUE", 
     ylab = "logFC", 
     col = "blue", 
     pch = 19)
grid()  # Adds grid lines to the plot
dev.off()  # Save and close the file


neuron_pseudo <-
  registration_pseudobulk(
    spe_ntc,
    var_registration = "neun_pos",
    var_sample_id = "sample_id",
    covars = c("age", "sex", "lot_num", "slide_id"),
    min_ncells = 10,
    pseudobulk_rds_file = here("processed-data", "image_processing", "enrichment", "neuron_pseudo.rds")
    )
	

dx_res <- registration_stats_enrichment(
  neuron_pseudo,
  block_cor = NaN,
  covars = c("age", "sex"),
  var_registration = "neun_pos",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

df <- read.csv(here("raw-data/images/SPG_Spot_Valid/NeuN/maddy.csv"), header = TRUE)

celltypes = c("Inhib", "Excit", "Excit_L3/4/5", "Excit_L3", "Excit_L4", "Excit_L6", "Excit_L5/6", "Excit_L5", "Excit_L2/3")
df_sub = df[df$cellType.target %in% celltypes,]
colnames(df_sub)[colnames(df_sub) == "gene"] <- "ensembl"

df_neuron  = merge(dx_res, df_sub, by = "ensembl")
Cneuron = cor(df_neuron$logFC_TRUE,df_neuron$logFC)

png(here("plots", "image_processing", "enrichment", "Neuron.png"), width = 800, height = 600)
plot(df_neuron$logFC_TRUE, df_neuron$logFC, 
     main = paste0("neuron ",dim(df_neuron)[[1]],",", Cneuron),
     xlab = "logFC_TRUE", 
     ylab = "logFC", 
     col = "blue", 
     pch = 19)
grid()  # Adds grid lines to the plot
dev.off() 