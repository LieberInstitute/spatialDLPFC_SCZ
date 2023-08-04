library(ggplot2)
library(PCAtools)
library(scran)
library(scater)

# TODO: 

# Read in SPE object (Optional) --------------------------------------------

method_batch_correct <- "harmony"

fld_batch_correct <- here("processed-data", "rds", "spe", "batch_corrected")

spe <- readRDS(
  here::here(fld_batch_correct,
             paste0("test_spe_", method_batch_correct,".rds")
  )
)



# Calculate Latent Dimensions ---------------------------------------------
set.seed(1)

## UMAP ------------------------------------------------------------------
# TODO: the dimred name can complicate
spe <- runUMAP(spe, dimred = "HARMONY", name = paste("UMAP.", method_batch_correct))
colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")




## PCA -------------------------------------------------------------------



## T-SNE  ----------------------------------------------------------------




# Qualitative Plots -------------------------------------------------------
## By Sample ID  -------------------------------------------------------
pdf(file = here::here("plots", "01_build_spe", "test_UMAP_harmony_sample_id.pdf"))
# Oh, my god, perfect match?
# TODO: check the presentation for a subset of samples that shares
#      different composition of morphology.
ggplot(
  data.frame(reducedDim(spe, "UMAP.HARMONY")),
  aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()
dev.off()
## By Seq Batch  -------------------------------------------------------

## By Slide ID  -------------------------------------------------------

