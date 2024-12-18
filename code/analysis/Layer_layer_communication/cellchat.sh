#!/bin/bash
#SBATCH --job-name=cellchat_layer_layer_comm
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=4            # Number of cores requested
#SBATCH --time=10:00:00 
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

## Load Modules
module load conda_R/4.4.x

## List current modules for reproducibility
module list

## Run code
Rscript cellchat.r


## Estimate JHPCE Memory
sstat -a -o MaxVMSize,AveVMSize,MaxRSS,AveRSS -j ${SLURM_JOB_ID}


echo "**** Job ends ****"
date
