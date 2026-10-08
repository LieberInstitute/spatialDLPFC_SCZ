#!/bin/bash
#SBATCH --job-name=08-cauchy_combination_per_sample
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH -n 1
#SBATCH --array=1-63%30
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x_%a.txt
#SBATCH --error=logs/%x_%a.txt    # file to collect standard output

# Per sample x trait: aggregate spot-level p-values within each spatial
# domain (spd_label) by Cauchy combination -> one p-value per domain.
# Requires: 07 for this sample.
# Output: {WORKDIR}/{sample}/cauchy_combination/{sample}_{trait}.Cauchy.csv.gz
#         (annotation, p_cauchy, p_median)

set -eo pipefail
mkdir -p logs
source gsMap_config.sh
log_job_info
activate_gsmap

sample_id=$(get_sample_id "${SLURM_ARRAY_TASK_ID}")
echo "Sample: ${sample_id}"

for trait in $(get_trait_names); do
  echo "---- trait ${trait} ----"
  gsmap run_cauchy_combination \
    --workdir "${WORKDIR}" \
    --sample_name "${sample_id}" \
    --trait_name "${trait}" \
    --annotation "${ANNOTATION}"
done

log_job_end
