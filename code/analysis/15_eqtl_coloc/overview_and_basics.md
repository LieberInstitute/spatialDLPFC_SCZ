## Project overview

We have a human brain (DLPFC) spatial transcriptomics/proteogenomics study with 63 donors where each donor has:

* gene expression measured on tissue sections (Visium HD-style data), producing a grid of spatial capture units (“spots” or “bins”) laid over the slice (SpatialExperiment  (SpE) R objects).

* spot labeling by 
  - "spatial domains" (SpD) 
  - spatial proteogenomics/microenvironment labels (SPG)

* Genotype data (genetic variants per donor)

The analytical goal is to perform eQTL mapping and colocalization analysis using the SpatialExperiment (SpE) gene expression: identify genetic variants associated with gene expression (eQTL), and then test whether those regulatory signals colocalize with GWAS loci (via coloc), in specific domains (SpD) and possibly microenvironment (SPG) defined contexts.


## Core spatial data terminology

  - Spots (spatial capture units): “spot” (or “bin” in HD) is a spatial measurement unit, not automatically a single cell. Each spot’s RNA can come from multiple nearby cells (2-10 cells). Spots live in physical 2D coordinate space of the tissue section/slice.

  - SpatialExperiment (SpE): The spatial data are stored as SpatialExperiment objects; it's like a SingleCellExperiment, but columns are spots/bins rather than cells, and it carries spatial coordinates + image-related metadata. Spot-level annotations (like layer or pathology labels) live as metadata columns attached to each spot.

  - SpD = spatial domains (anatomical layer / lamina) ; SpD labels assign each spot to an anatomical domain, most commonly the cortical layers. Interpretation: *"Where is this spot in the layered structure of cortex?"*
  In this study we have: SpD01_WMtz, SpD02_L3-4, SpD03_L6, SpD04_WM, SpD05_L5, SpD06_L2-3, SpD07_L1-M where WMtz=White_Matter_transition_zone, WM=White_Matter, L=Layer, M=Meninges

  - SPG = spatial proteogenomic (microenvironment/IF labeling); SPG labels assign each spot to a  group, derived from tissue imaging signals (IF) and spatial proximity rules. PNNs (perineuronal nets) are a biological feature that can be measured by staining/imaging (e.g., a PNN marker channel). The imaging signal is aligned to the spatial grid, and spots are classified based on PNN presence/intensity and neighborhood. In this study, the microenvironment labels are: neuropil, PNN, neuronal, vasculature.


## Using tensorQTL in this project

Genotype varies by donor, so eQTL testing must be based on one expression value per donor per gene (for the phenotype used in a given run). TensorQTL cannot make use of multiple expression data samples for the same genotype, so expression data aggregation is necessary. We run separate eQTL mapping analyses per "context" (e.g., per SpD pseudobulk data), using one sample per donor per run.

* SpD pseudobulk (domain-specific expression): ignoring SPG labels, run eQTL separately for each SpD context
* SPG pseudobulk (ignoring layers): aggregate expression data by microenvironment, ignoring SpD assignment
* SpD x SPG (layer-by-microenvironment): for each donor and each (layer, microenvironment) combination: aggregate spots in that intersection. This can explode the number of contexts and destroy effective sample size; usually only feasible for a few well-populated combinations.
