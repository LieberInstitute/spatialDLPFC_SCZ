#!/bin/bash
#SBATCH --job-name=plot_outlier_spots
#SBATCH --mem=40G
#SBATCH -n 1
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
Rscript plot_outlier_spots.R


## Memeory stat
sstat -a -o JobID,MaxVMSizeNode,MaxVMSize,AveVMSize,MaxRSS,AveRS S,MaxDiskRead,MaxDiskWrite,AveCPUFreq,TRESUsageInMax -j ${SLURM_JOB_ID}

echo "**** Job ends ****"
date

