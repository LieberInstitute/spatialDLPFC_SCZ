#!/bin/bash
#SBATCH --job-name=PRECAST_test
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=4            # Number of cores requested
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x.txt
#SBATCH --error=logs/%x.txt    # file to collect standard output


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
Rscript PRECAST.R


## Estimate JHPCE Memory
sstat -o MaxVMSize,AveVMSize,MaxRSS,AveRSS -j ${SLURM_JOB_ID}


echo "**** Job ends ****"
date
