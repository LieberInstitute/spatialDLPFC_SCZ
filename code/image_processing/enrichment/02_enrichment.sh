#!/bin/bash
#SBATCH --mem=400G
#SBATCH --job-name=enrichmentStats
#SBATCH -o logs/enrichmentStats.txt 
#SBATCH -e logs/enrichmentStats.txt 

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"s
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load MATLAB
module load conda_R/devel

Rscript 02_enrichment.R


echo "**** Job ends ****"
date