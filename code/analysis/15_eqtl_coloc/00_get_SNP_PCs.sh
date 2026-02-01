#!/bin/bash
prefix=genotypes/plink2/merged_maf05
## check if genotypes/$prefix.pgen etc exist, error if not
if [ ! -f ${prefix}.pgen ]; then
  echo "ERROR: ${prefix}.pgen not found, please check the prefix"
  exit 1
fi
## First use indep-pairwise to make sure that PCs are not biased by
## linkage disequilibrium (LD)
## (assumes basic qc has been performed)
plink2 --pfile $prefix --indep-pairwise 50 5 0.2 --out ${prefix}_ldpruned
## --indep-pairwise 50 5 0.2: performs LD pruning with a window size of 50 SNPs,
##                   shifting the window by 5 SNPs at a time, and removing SNPs with an r² greater than 0.2.

## After pruning, use the pruned SNPs to compute PCs:
plink2 --pfile $prefix --extract ${prefix}_ldpruned.prune.in \
       --pca 5  --out ${prefix}_pca
## ${prefix}_pca.eigenvec file has the PCs, the first column is the sample ID
# #IID	PC1	PC2	PC3	PC4	PC5
# Br1113	0.0339709	0.196971	-0.13841	0.312701	-0.0992901
# Br1753	0.00134058	-0.0269248	0.0823359	-0.119852	-0.104399
