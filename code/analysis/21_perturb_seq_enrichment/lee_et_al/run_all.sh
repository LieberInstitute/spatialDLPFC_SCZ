#!/bin/bash
# Run the Lee et al. (2026) perturbation enrichment analysis in order.
# Usage: bash run_all.sh   (from this folder; scripts resolve paths with here())
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs

for script in \
  00-prep_lee_gene_sets.R \
  01-prep_our_rankings.R \
  02-fgsea_lee_sets_vs_our_rankings.R \
  03-plot_fgsea_lee_sets_vs_our_rankings.R \
  04-fgsea_our_sets_vs_lee_rankings.R \
  05-plot_fgsea_our_sets_vs_lee_rankings.R \
  06-plot_ortholog_matching.R \
  07-fig6C_check_all_DEGs.R \
  08-fgsea_lee_pooled_sets_vs_our_rankings.R \
  09-plot_fgsea_lee_pooled_sets_vs_our_rankings.R \
  10-fgsea_our_pooled_sets_vs_lee_rankings.R \
  11-plot_fgsea_our_pooled_sets_vs_lee_rankings.R; do
  echo "**** Running ${script} ****"
  date
  Rscript "${script}" > "logs/${script%.R}.txt" 2>&1
done

# Reverse pooled analysis with our DEGs at nominal p < 0.05
echo "**** Running 10/11 with nominal p < 0.05 ****"
date
Rscript 10-fgsea_our_pooled_sets_vs_lee_rankings.R p05 > logs/10-fgsea_our_pooled_sets_vs_lee_rankings_nomp05.txt 2>&1
Rscript 11-plot_fgsea_our_pooled_sets_vs_lee_rankings.R p05 > logs/11-plot_fgsea_our_pooled_sets_vs_lee_rankings_nomp05.txt 2>&1

echo "**** All steps finished ****"
date
