library(batchelor)
library(scuttle)

# Create a list of sample-specific SPEs
spe_lst <- map()

# Feed the list to multiBatchNorm
multiBatchNorm()


# Create modelGeneVar for each item of the list

HVG_lst <- spe_lst |> map( .f = modelGeneVar) |> 
  #TODO: flatten

# Use combineVar to combine these variables
combined.dec <- combineVar(HVG_lst)
chosen.hvgs <- getTopHVGs(combined.dec, n=5000)

# Validation of batch effect
combined <- correctExperiments(spe_lst,
                               PARAM=NoCorrectParam()) # Note not corrected 

library(scater)
set.seed(100)
combined <- runPCA(combined, subset_row=chosen.hvgs)
combined <- runTSNE(combined, dimred="PCA")
plotTSNE(combined, colour_by="batch")


# Fast MNN
set.seed(101)
f.out2 <- fastMNN(combined, batch=combined$batch, subset.row=chosen.hvgs)


set.seed(103)
f.out <- runTSNE(f.out, dimred="corrected")
plotTSNE(f.out, colour_by="batch")