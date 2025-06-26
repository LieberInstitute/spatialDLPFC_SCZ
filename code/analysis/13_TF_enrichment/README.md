We manually ran ChEA3 (Version 3) with our layer-restricted DEGs (nomical p_value < 0.05). 

Steps:
1. run 00-parepare_files_for_ChEA3.R to create up-/down-regulated DEGs per spatial domains.
2. Manually run each DEGs per domains using https://maayanlab.cloud/chea3/ portal.    
  - Retrospectively thinking, this could be done with code if the interest is only to get the meanRank results. See https://maayanlab.cloud/chea3/index.html#content4-z 
3. run processing scripts after the TF results. 