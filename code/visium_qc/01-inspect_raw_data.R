# vis_grid_gene(
#   spe = spe,
#   geneid = "Ndata",
#   pdf = here::here("plots", "01_build_spe", "seg_count.pdf"),
#   assayname = "counts",
#   auto_crop = FALSE,
#   spatial = FALSE
# )

# mean(spe$count)
# pdf(here::here("plots", "01_build_spe", "cells_per_spot.pdf"))
# boxplot(spe$count)
# dev.off()

## Remove genes with no data
# no_expr <- which(rowSums(counts(spe)) == 0)

# length(no_expr)
# [1] 6936
# length(no_expr) / nrow(spe) * 100
# [1] 18.9503



