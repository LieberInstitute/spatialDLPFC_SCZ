#!/bin/bash
#SBATCH --job-name=02-project_SNAP_scores
#SBATCH --mem=40G
#SBATCH --time=00:45:00
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x.txt
#SBATCH --error=logs/%x.txt    # file to collect standard output

mkdir -p logs

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_CLUSTER_NAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load Modules
module load conda_R/devel

## List current modules for reproducibility
module list

## Run code
Rscript 02-project_SNAP_scores.R

echo "**** Job ends ****"
date
