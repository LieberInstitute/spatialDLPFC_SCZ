
# Load Library ------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(sessioninfo)
  library(tidyverse)
  library(spatialLIBD)
  library(tidyverse)
})


# Load Data -----
## Load NMF Res --------
model <- readRDS(
  here(
    "processed-data/rds/NMF/",
    "test_NMF_all_k50.rds"
  )
)

## Load Spe ---------------
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster")

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform",".rds")
)

spe <- readRDS(
  path_PRECAST_int_spe
)



# For test only
# spe_backup <- spe
# spe <- spe_backup

spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )]

# k <- 2
k <- 50
patterns <- t(model$h) # these are the factors
colnames(patterns) <- paste0("NMF", 1:k)
colData(spe) <- cbind(colData(spe), patterns)


# Spot Plot of NMF Patterns ----
.sample_ordered <- metadata(spe)$dx_df |> arrange(dx, sample_id) |> pull(sample_id)


for(.fac in paste0("NMF", 1:k)){
  # .fac <- "NMF1"
  vis_grid_gene(
    spe,
    geneid = .fac,
    pdf_file = here(paste0("plots/NMF/test_spot_plot_", .fac, ".pdf")),
    sample_order = .sample_ordered,
    point_size = 0.8,
    spatial = FALSE
  )
}


# Correlation of Factors ----

mat_NMF_cor <- cor(patterns)
pdf(here("plots/NMF/NMF_corr_mat.pdf"), width = 10, height = 10)
heatmap(mat_NMF_cor, symm = TRUE, scale = "none")
dev.off()

pdf(here("plots/NMF/abs_NMF_corr_mat.pdf"), width = 15, height = 15)
heatmap(abs(mat_NMF_cor), symm = TRUE, scale = "none")
dev.off()






# library(tidyverse)
# subset_id <- metadata(spe)$dx_df |> group_by(dx) |>
#   slice_head(n=2) |> ungroup() |>
#   pull(sample_id)
# 
# spe <- spe[ ,spe$sample_id %in% subset_id]

# ggplot() +
#   geom_violin(aes(y = spe$NMF1))

# spe$NMF1 |> density() |> plot()


# library(spatialLIBD)
# vis_grid_gene(
#   spe, geneid = "NMF1", spatial = FALSE,
#   pdf_file = here(
#     "plots/NMF",
#     paste0("test_NMF_k",k,"_F1", ".pdf")
#   )
# )
# 
# 
# 
# vis_grid_gene(
#   spe, geneid = "NMF2", spatial = FALSE,
#   pdf_file = here(
#     "plots/NMF",
#     paste0("test_NMF_k",k,"_F2", ".pdf")
#   )
# )
# It feel like the NMF are very zero inflated.
density(spe$NMF3) |> plot()


# TODO: visualize what the two factors are 
## Scatter plot....

# t.test(spe$NMF1~spe$dx)
# t.test(spe$NMF2~spe$dx)
pdf(here("plots/NMF/test_nmf_k50_all_spots_density_comp.pdf"))
paste0("NMF", 1:k) |> 
  walk(.f = function(nmf_name){
    (ggplot() +
      geom_density(aes(x = spe[[nmf_name]], group = spe$dx)) +
       labs(title = nmf_name)) |> 
      print()
    
  })
dev.off()

col_df <- colData(spe) |> data.frame()
all_spots_test <- paste0("NMF", 1:k) |> 
  map(
    ~t.test(col_df[[.x]] ~ col_df$dx)
  )

lapply(all_spots_test, 
      FUN = function(.res){
        # browser()
        .res$p.value
      }
) |> unlist() |> 
  set_names(paste0("NMF", 1:k))

pdf(here("plots/NMF/test_nmf_k50_all_spots_registration.pdf"))

for(spd_var in paste0("PRECAST_", 2:9)){
  # browser()
  # spd_var <- "PRECAST_9"
  # stopifnot(!all(is.numeric(spe[[spd_var]])))
  dn_mat_spd <- model.matrix(~factor(spe[[spd_var]])-1) 
  # kronecker(dn_mat_spd)
  
  # dn_mat_spd |> str()
  
  
  
  # ggplot() +
  #   geom_density(aes(x = spe$NMF2, group = spe$dx))
  
  reg_mat <- cor(col_df |> select(starts_with("NMF")), dn_mat_spd)
  
  # reg_mat_long <- reg_mat |>   reshape2::melt()
  
  # ggplot(data = reg_mat_long) +
  #   geom_tile(aes(Var1, Var2, fill = value)) +
  #   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  #   theme_minimal() +
  #   labs(title = "NMF Registering to PRECAST SpD9",
  #        x = "NMF",
  #        y = "SpD")
  # 
  
  
  
  heatmap(reg_mat,#Rowv = NA, Colv = NA,
          scale = "none",
          # col = colorRampPalette(c("blue", "white", "red"))(100),
          # symm = TRUE,  # Show the matrix symmetrically
          main = paste0("NMF Registering to ", spd_var),
          xlab = "SpD",
          ylab = "NMF"
  ) |> print()
  
}




# NMF registration to sample ----------------------------------------------
dn_mat_smp <- model.matrix(~spe$sample_id-1)
reg_mat <- cor(col_df |> select(starts_with("NMF")), dn_mat_smp)
heatmap(reg_mat,Colv = NA, #Rowv = NA, 
        scale = "none",
        # col = colorRampPalette(c("blue", "white", "red"))(100),
        # symm = TRUE,  # Show the matrix symmetrically
        main = "NMF Registration to sample",
        xlab = "Sample",
        ylab = "NMF"
) |> print()

dev.off()






# TODO: Test if the factor is associated with 

print(
  "Finish the analysis"
)
# Session Info ------------------------------------------------------------
sessioninfo::session_info()

