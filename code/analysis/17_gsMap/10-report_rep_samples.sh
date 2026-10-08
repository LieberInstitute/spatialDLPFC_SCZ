#!/bin/bash
#SBATCH --job-name=10-report_rep_samples
#SBATCH --mem=40G
#SBATCH --time=4:00:00
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x.txt
#SBATCH --error=logs/%x.txt    # file to collect standard output

# Optional: gsMap's interactive HTML report (spot-level -log10 p map,
# domain-level Cauchy p, top GSS genes correlated with the trait signal)
# for the representative samples, as a quick QC/exploration of each trait.
# Requires: 08 for these samples.
# Output: {WORKDIR}/{sample}/report/{trait}/

set -eo pipefail
mkdir -p logs
source gsMap_config.sh
log_job_info
activate_gsmap

for sample_id in "${REP_SAMPLES[@]}"; do
  for trait in $(get_trait_names); do
    echo "---- sample ${sample_id}, trait ${trait} ----"
    gsmap run_report \
      --workdir "${WORKDIR}" \
      --sample_name "${sample_id}" \
      --trait_name "${trait}" \
      --annotation "${ANNOTATION}" \
      --sumstats_file "$(get_sumstats_file "${trait}")" \
      --top_corr_genes 50
  done
done

log_job_end
