#!/bin/bash
#SBATCH --job-name=01-export_spe_to_h5ad
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x.txt
#SBATCH --error=logs/%x.txt    # file to collect standard output

mkdir -p logs
source gsMap_config.sh
log_job_info

## Load Modules
module load conda_R/devel

## List current modules for reproducibility
module list

## Run code
Rscript 01-export_spe_to_h5ad.R

log_job_end
