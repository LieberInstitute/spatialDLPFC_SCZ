
  
# Run the Python script
python PRScs/PRScs.py --ref_dir=PRScs/Reference/ldblk_1kg_eur \
--bim_prefix=/data/vbalakundi/Spatial_DLPFC/data/genotype/EUR/EUR.merged_R.8_MAF.01.RSann.noDups.QC.01 \
--sst_file=/data/vbalakundi/GWAS/SCZ/Trubetskoy_2022_Nature/PGC3_SCZ_wave3.european.autosome.public.v3_PRScs_masked_chr3_chr6.txt \
--n_gwas=127906 \
--out_dir=/data/vbalakundi/Spatial_DLPFC/Results/PRScs/SCZ_PRScs_0.01_masked_chr3_chr6/EachChromosome/



# Combine per-chromosome scores for Masked chr3 + chr6

cat "${EUR_PRScs}/SCZ_PRScs_0.01_masked_chr3_chr6/EachChromosome/"*.txt > "${EUR_PRScs}/SCZ_PRScs_0.01_masked_chr3_chr6/SCZ_PRScs_0.01_masked_chr3_chr6_Scores.txt"
 

# e.r. Run PRS scoring for Masked chr3 + chr6

plink \

  --bfile "${EUR_Genotype}/EUR.merged_R.8_MAF.01.RSann.noDups.QC.01" \

  --score "${EUR_PRScs}/SCZ_PRScs_0.01_masked_chr3_chr6/SCZ_PRScs_0.01_masked_chr3_chr6_Scores.txt" 2 4 6 sum \

  --out "${EUR_PRScs}/SCZ_PRScs_0.01_masked_chr3_chr6/SCZ_PRScs_0.01_masked_chr3_chr6_EUR_Scores"
 


 