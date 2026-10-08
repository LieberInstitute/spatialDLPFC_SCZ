#!/bin/bash
#SBATCH --job-name=09-cauchy_combination_across_samples
#SBATCH --mem=40G
#SBATCH --time=2:00:00
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x.txt
#SBATCH --error=logs/%x.txt    # file to collect standard output

# Per trait: pool spot-level p-values for each spatial domain across samples
# by Cauchy combination, for (a) all 63 samples, (b) NTC samples only and
# (c) SCZ samples only.
# Requires: 07 for all samples.
# Output: {CAUCHY_ACROSS_DIR}/{trait}_{group}.Cauchy.csv.gz
#         (annotation, p_cauchy, p_median)

set -eo pipefail
mkdir -p logs
source gsMap_config.sh
log_job_info
activate_gsmap

mkdir -p "${CAUCHY_ACROSS_DIR}"

declare -A group_list=(
  [all]="${SAMPLE_LIST}"
  [ntc]="${SAMPLE_LIST_NTC}"
  [scz]="${SAMPLE_LIST_SCZ}"
)

for trait in $(get_trait_names); do
  for group in all ntc scz; do
    echo "---- trait ${trait}, group ${group} ----"
    mapfile -t sample_names < "${group_list[${group}]}"
    echo "N samples: ${#sample_names[@]}"

    gsmap run_cauchy_combination \
      --workdir "${WORKDIR}" \
      --sample_name_list "${sample_names[@]}" \
      --trait_name "${trait}" \
      --annotation "${ANNOTATION}" \
      --output_file "${CAUCHY_ACROSS_DIR}/${trait}_${group}.Cauchy.csv.gz"
  done
done

log_job_end
