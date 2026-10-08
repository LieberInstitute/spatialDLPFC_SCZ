#!/bin/bash
#SBATCH --job-name=03-create_slice_mean
#SBATCH --mem=80G
#SBATCH --time=6:00:00
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x.txt
#SBATCH --error=logs/%x.txt    # file to collect standard output

# Compute the across-slice geometric mean of gene expression ranks over all
# 63 samples. Passed to run_latent_to_gene (05) via --gM_slices, so gene
# specificity scores (GSS) are on a common scale across samples and are
# comparable between donors/diagnoses.

set -eo pipefail
mkdir -p logs
source gsMap_config.sh
log_job_info
activate_gsmap

mkdir -p "$(dirname "${SLICE_MEAN_FILE}")"

mapfile -t sample_names < "${SAMPLE_LIST}"
mapfile -t h5ad_files < "${H5AD_DIR}/h5ad_list.txt"

# error prevention: lists are aligned
if [ "${#sample_names[@]}" -ne "${#h5ad_files[@]}" ]; then
  echo "sample_list.txt and h5ad_list.txt differ in length" >&2
  exit 1
fi

gsmap create_slice_mean \
  --sample_name_list "${sample_names[@]}" \
  --h5ad_list "${h5ad_files[@]}" \
  --slice_mean_output_file "${SLICE_MEAN_FILE}" \
  --data_layer "${DATA_LAYER}"

log_job_end
