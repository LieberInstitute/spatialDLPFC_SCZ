# spatialDLPFC_SCZ



## Goal: Mapping spatially organized molecular and genetic signatures of schizophrenia across multiple scales in the human prefrontal cortex.

![Experimental Overview](https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/gh-pages/img/overview.png?raw=true)  

## Study Design  

Molecular mechanisms underlying dorsolateral prefrontal cortex (dlPFC) dysfunction in schizophrenia (SCZ) are poorly understood. dlPFC cell types are spatially organized across the six layers into functional microcircuits, which regulate cognitive and emotional processes that are implicated in SCZ. While regional specificity across cortical layers and cell types has been demonstrated for some SCZ-linked genes, spatially-resolved transcriptomics (SRT) can more definitively map molecular associations of disease. We investigated spatial gene expression changes in the human dlPFC from neurotypical control (n=31) and SCZ (n=32) brain donors using the Visium platform with incorporation of immunostaining to label perineuronal nets, neurons and vasculature. SCZ-associated DEGs were then mapped across these labeled cellular microenvironments and cortical layers. Major transcriptional alterations were identified in synaptic and neuroimmune pathways, which we localized to neuropil and glia-enriched domains. Integrative analysis with bulk and single-cell RNA studies highlighted distinct roles for spatially-localized glial cell populations, and identified enrichment of novel DEGs in endothelial cells and microglia. These findings were supported by enrichment of SCZ genetic risk across similar domains, and association of laminae-specific transcription factors to SCZ risk variants. Cellular resolution SRT using Xenium extended laminar disease associations to spatially-organized cell types. The findings highlight unique advantages of SRT to identify novel SCZ-related biology, particularly for cellular populations and compartments that may be missed with alternative technologies. We provide data resources to enable sample-level visualization of the spatial transcriptomics data (Samui) and to explore SCZ-associated DEGs across spatial domains (iSEE-layer adjusted) and cellular microenvironments (iSEE-Neuropil, iSEE-Neuron, iSEE-vasculature, iSEE-PNNs).  

For the associated Xenium repository for this project, please visit [spatialDLPFC_SCZ_XENIUM](https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM).

## Interactive Websites

All of these interactive websites are powered by open source software,
namely:  

- 🔍 [`samui`](http://dx.doi.org/10.1017/S2633903X2300017X)  
- 👀 [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)  

We provide the following interactive websites, organized by dataset with
software labeled by emojis:  

- Visium (n = 63)
  - 🔍 [DLPFC SCZ Samui browser - NTC, F (N=15)](https://samuibrowser.com/from?url=data.libd.org/samuibrowser/&s=Br2720_V13M06-280_C1&s=Br5367_V12F14-053_D1&s=Br5400_V13M06-281_A1&s=Br6297_V13M06-282_C1&s=Br2719_V12F14-053_A1&s=Br1113_V13M06-279_A1&s=Br8042_V13M06-279_B1&s=Br5314_V13F27-294_C1&s=Br8325_V13F27-295_B1&s=Br8492_V13M06-343_B1&s=Br5622_V13F27-293_A1&s=Br1753_V13F27-296_C1&s=Br5931_V13M06-340_A1&s=Br5228_V13M06-344_D1&s=Br8667_V13M06-342_D1)
  - 🔍 [DLPFC SCZ Samui browser - NTC, M (N=16)](https://samuibrowser.com/from?url=data.libd.org/samuibrowser/&s=Br5182_V12F14-053_C1&s=Br1092_V13F27-295_A1&s=Br1204_V13F27-336_B1&s=Br5395_V13M06-281_B1&s=Br5472_V12D07-334_B1&s=Br6432_V12D07-334_D1&s=Br5276_V13M06-343_A1&s=Br5436_V13M06-282_D1&s=Br6257_V13M06-280_D1&s=Br8433_V13M06-342_C1&s=Br5599_V13F27-293_B1&s=Br5639_V13F27-294_D1&s=Br6389_V13F27-296_D1&s=Br1214_V13M06-344_C1&s=Br8406_V13M06-340_B1&s=Br8218_V13F27-336_A1)
  - 🔍 [DLPFC SCZ Samui browser - SCZ, F (N=14)](https://samuibrowser.com/from?url=data.libd.org/samuibrowser/&s=Br1958_V12F14-057_B1&s=Br5590_V13M06-280_A1&s=Br8537_V13M06-279_C1&s=Br8514_V13M06-279_D1&s=Br8772_V13F27-295_D1&s=Br2242_V13M06-281_C1&s=Br5746_V13F27-294_A1&s=Br6560_V13M06-282_A1&s=Br5588_V13F27-293_C1&s=Br5973_V13M06-343_D1&s=Br8207_V13M06-342_B1&s=Br1139_V13F27-296_A1&s=Br2470_V13M06-340_D1&s=Br6369_V13M06-344_B1)
  - 🔍 [DLPFC SCZ Samui browser - SCZ, M  (N=18)](https://samuibrowser.com/from?url=data.libd.org/samuibrowser/&s=Br1526_V12F14-057_A1&s=Br2039_V12F14-057_C1&s=Br2378_V12F14-057_D1&s=Br5373_V12D07-334_A1&s=Br5636_V13M06-280_B1&s=Br5789_V12D07-334_C1&s=Br5888_V13M06-281_D1&s=Br5928_V13F27-294_B1&s=Br1682_V13M06-340_C1&s=Br5224_V13M06-342_A1&s=Br5573_V13M06-343_C1&s=Br6032_V13F27-293_D1&s=Br6437_V13M06-282_B1&s=Br1923_V13F27-336_C1&s=Br2627_V13F27-295_C1&s=Br1578_V13F27-296_B1&s=Br2421_V13F27-336_D1&s=Br6496_V13M06-344_A1)
    -  Provides interactive spot-level visualization of Visium data.
- Pseudobulk and micro environment (n = 63)
  - 👀 [pseudobulked iSEE app](https://interactive.libd.org/spatialDLPFC_scz-spatial_domains/)
  - 👀 [Neuronal iSEE app](https://interactive.libd.org/spatialDLPFC_SCZ-neuronal/)
  - 👀 [Neuropil iSEE app](https://interactive.libd.org/spatialDLPFC_SCZ-neuropil/)
  - 👀 [PNN iSEE app](https://interactive.libd.org/spatialDLPFC_scz-PNN/)
  - 👀 [Vasc iSEE app](https://interactive.libd.org/spatialDLPFC_scz-vasc/)


## Data Access
Public [globus endpoint](https://research.libd.org/globus/) to access R objects associate with apps for this project.  
All data, including raw FASTQ files and SpaceRanger processed data outputs, can be accessed via Gene Expression Omnibus (GEO) under accessions [GSE307403](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307403) and [GSE307404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307404)
Zenodo Archive for this project can be found at [10.5281/zenodo.18663825](https://doi.org/10.5281/zenodo.18663825).

## Background:  
[https://www.sciencedirect.com/science/article/pii/S092099641500002X?via%3Dihub](https://www.sciencedirect.com/science/article/pii/S092099641500002X?via%3Dihub)  
[https://www.frontiersin.org/articles/10.3389/fnmol.2018.00270/full](https://www.frontiersin.org/articles/10.3389/fnmol.2018.00270/full)  
[https://www.nature.com/articles/s41583-019-0196-3](https://www.nature.com/articles/s41583-019-0196-3)  
[https://pubmed.ncbi.nlm.nih.gov/34072323/](https://pubmed.ncbi.nlm.nih.gov/34072323/)  

## Contact

We value public questions, as they allow other users to learn from the
answers. If you have any questions, please ask them at
[LieberInstitute/spatialDLPFC_SCZ/issues](https://github.com/LieberInstitute/spatialDLPFC_SCZ/issues)
and refrain from emailing us. Thank you again for your interest in our
work!
