#!/bin/bash
#SBATCH --job-name=01-split_h5ad_by_sample
#SBATCH --mem=40G
#SBATCH --time=2:00:00
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=boyi.guo@utah.edu
#SBATCH --output=logs/%x.txt
#SBATCH --error=logs/%x.txt    # file to collect standard output

# Split processed-data/samui/spg_kept_only.h5ad (full dataset, from
# code/samui/spe_to_anndata.R) into one gsMap-ready h5ad per sample.
# Requires the gsMap conda env from 00-setup_gsMap_env.sh.

set -eo pipefail
mkdir -p logs
source gsMap_config.sh
log_job_info
activate_gsmap

python 01-split_h5ad_by_sample.py "${SOURCE_H5AD}" "${H5AD_DIR}"

log_job_end
