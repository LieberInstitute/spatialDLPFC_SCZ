## Files:

- cellchat.r: running comparative cell-cell communication between the diagnosis groups. It first run the cell-chat pipeline on the subset of samples from each diagnosis separately (stored), and then use xxx fucntion to merge the two cellchat object and identify the differential pattern and measure the strength of the evidence.

- DEG_only: Folder containing code to characterize dx-DEGs in the data-driven layer-layer communication analysis.
  - [DEG_only/understanding_pathways.R](DEG_only/understanding_pathways.R) contain characterization and plotting. 
  - *.csv containing the pathways/ligand receptor pairs where DEGs serves as ligand or receptor genes.

- `/Users/bguo6/GitHub/spatialDLPFC_SCZ/code/analysis/Layer_layer_communication/vulnerablility_charact` studies the vulnerable layers, especially L2/3 and L6/WM differential pathways.

- Validation
  - cellchat_validate_NTC.qmd/.pdf: Visualizes the results generated from `cellchat.r`, provides summary report and 

  As a part of validation of the concept of Layer-layer communication analysis as well as the sensitivity/specificity of cellchat.r, we run some validation analysis. 

  - cell_chat_validate_NTC.qmd/.pdf: a summary report focuing on the communication pattern among only NTC samples. It leverages the results from `cellchat.r`


  - call_martinotti-in_L5.r is an exploratory analysis to see if we use some martinotti marker gene from mouse model to call spots in L5. some **Viz** is provided too. 

  - Ligand_expression_in_L5.qmd/.pdf: a more complete report explaining and docmument the exploratory analysis of martinotti+ spot calling.

  - the cellchat_martinotti.r + cellchat_validate_NTC_martinotti.qmd/.pdf: re-run the analysis and report among L5 (**both martinotti+ and martinotti - spots**) + L1 spots only. 

  - cellchat_martinotti_only.r + cellchat_validate_NTC_martinotti_only.qmd/.pdf:  re-run the analysis and report among a subset of L5 spots (**martinotti+ spots only**) + L1 spots only. 

- Archived
  - liana_draft.r: Layer-to-Layer communication analysis using `liana`. The analysis was run by each diagnosis group respectively under the variable `.dx`
  - new_compare: a folder containing codes which *runs cellchat on each sample*. *Sample-specific pathway information is used to do a t-test* to see if there's any pathway are different between the two diagnosis. So far, the t-test only yields negative results. and hence, deprecate this effort.