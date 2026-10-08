#!/bin/bash
#SBATCH --job-name=04-find_latent_representations
#SBATCH --mem=30G
#SBATCH --time=4:00:00
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --array=1-63%20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x_%a.txt
#SBATCH --error=logs/%x_%a.txt    # file to collect standard output

# Per sample: train gsMap's graph neural network (GNN-VAE) on expression +
# spatial neighbourhood, with spd_label as the supervising annotation.
# Output: {WORKDIR}/{sample}/find_latent_representations/{sample}_add_latent.h5ad
# (obsm['latent_GVAE']).

set -eo pipefail
mkdir -p logs
source gsMap_config.sh
log_job_info
activate_gsmap

sample_id=$(get_sample_id "${SLURM_ARRAY_TASK_ID}")
echo "Sample: ${sample_id}"

gsmap run_find_latent_representations \
  --workdir "${WORKDIR}" \
  --sample_name "${sample_id}" \
  --input_hdf5_path "${H5AD_DIR}/${sample_id}.h5ad" \
  --annotation "${ANNOTATION}" \
  --data_layer "${DATA_LAYER}"

log_job_end
