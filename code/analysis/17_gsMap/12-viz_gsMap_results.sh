#!/bin/bash
#SBATCH --job-name=12-viz_gsMap_results
#SBATCH --mem=60G
#SBATCH --time=2:00:00
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
Rscript 12-viz_gsMap_results.R

log_job_end
