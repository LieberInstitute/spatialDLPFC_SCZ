# NOTE: Code adapted from https://github.com/dssikdar/C2Cv0


# Load Package ----
# library(Seurat)
# library(SeuratDisk)
# library(NMF)
library(here)
library(ggalluvial)
library(ComplexHeatmap)
library(CellChat)
library(patchwork)
library(graphics)
# library(Matrix)
library(circlize)
# library(colorspace)
# library(pracma)
# library(glue)

# options(stringsAsFactors = FALSE, repr.plot.width = 12, repr.plot.height = 9, repr.plot.res = 300)

# Load helper function.
source(here(
  "code/analysis/SPG_analysis/ccc_neun_claudin",
  "3_A-Sender_Receiver_Chord.r"
))

# Load Data ---
## Load merged object ----
cellchat <- readRDS(
  here(
    "processed-data/rds/spg_ccc",
    "spd_neun_ex_merged_cellchat.rds"
  )
)

## Load NTC and SCZ cellchat objects ----
ntc_cellchat <- readRDS(here(
  "processed-data/rds/spg_ccc",
  "spd_neun_ex_ntc_cellchat.rds"
))
scz_cellchat <- readRDS(
  here(
    "processed-data/rds/spg_ccc",
    "spd_neun_ex_scz_cellchat.rds"
  )
)

object.list <- list(ntc = ntc_cellchat, scz = scz_cellchat)

## Create merged merged dataset ---
measure <- "weight"

### (Deprecated) NOte this section doesn't seem to be that useful ---
# obj1 <- object.list[[1]]@net[[measure]]              # NTC
# obj2 <- object.list[[2]]@net[[measure]]              # SCZ
# Turn the matrix to long form
# melt_1 <- reshape2::melt(obj1, value.name="count")
# melt_2 <- reshape2::melt(obj2, value.name="count")
# sum1 <- sum(melt_1$count)
# sum2 <- sum(melt_2$count)

# # balance
# normalization_NTCvSCZ = sum1 / sum2
# # normalization = sum(object.list[[1]]@netP$prob) / sum(object.list[[2]]@netP$prob)
# net.diff_NTCvSCZ <- obj2 * (normalization_NTCvSCZ) - obj1

# Make plots
## Idenitify all pathways ----
all_pathways <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE, return.data = TRUE)[["signaling.contribution"]]


all_pw_names <- all_pathways$name |> unique()

selected_pw_names <- c("CCK", "GRN", "2-AG", "Glutamate", "GABA-B", "GABA-A", "NRXN")
stopifnot(all(selected_pw_names %in% all_pw_names))

## Make plot -----
# This is not the plot.
# gg2 <- my_netVisual_aggregate2(net.diff_NTCvSCZ, layout="chord",  signaling.name = all_pw_names[1])

# Red means enriched SCZ

pdf(here("code/analysis/SPG_analysis/ccc_ex_neun/chord_plot_diff_path.pdf"))
for (.path in selected_pw_names) {
  print(my_netVisual_aggregate(ntc_cellchat, scz_cellchat, signaling = .path, measure = "weight", layout = "chord"))
}
dev.off()
