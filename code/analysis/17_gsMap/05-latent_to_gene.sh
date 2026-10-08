#!/bin/bash
#SBATCH --job-name=05-latent_to_gene
#SBATCH --mem=40G
#SBATCH --time=4:00:00
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --array=1-63%20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x_%a.txt
#SBATCH --error=logs/%x_%a.txt    # file to collect standard output

# Per sample: compute per-spot gene specificity scores (GSS) from the
# latent-space + spatial neighbourhood of each spot, using the across-sample
# slice mean (03) so GSS are comparable across samples.
# Requires: 03 (slice mean) and 04 (latent representation) for this sample.
# Output: {WORKDIR}/{sample}/latent_to_gene/{sample}_gene_marker_score.feather

set -eo pipefail
mkdir -p logs
source gsMap_config.sh
log_job_info
activate_gsmap

sample_id=$(get_sample_id "${SLURM_ARRAY_TASK_ID}")
echo "Sample: ${sample_id}"

# error prevention
test -s "${SLICE_MEAN_FILE}"

gsmap run_latent_to_gene \
  --workdir "${WORKDIR}" \
  --sample_name "${sample_id}" \
  --annotation "${ANNOTATION}" \
  --latent_representation latent_GVAE \
  --num_neighbour "${NUM_NEIGHBOUR}" \
  --num_neighbour_spatial "${NUM_NEIGHBOUR_SPATIAL}" \
  --gM_slices "${SLICE_MEAN_FILE}"

log_job_end
