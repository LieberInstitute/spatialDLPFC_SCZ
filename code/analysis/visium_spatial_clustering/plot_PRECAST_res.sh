#!/bin/bash
#SBATCH --job-name=plot_PRECAST_res
#SBATCH --mem=30G
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1            # Number of cores requested
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
module load conda_R/devel

## List current modules for reproducibility
module list

## Run code
Rscript plot_PRECAST_res.R


## Estimate JHPCE Memory
sstat -a -o MaxVMSize,AveVMSize,MaxRSS,AveRSS -j ${SLURM_JOB_ID}


echo "**** Job ends ****"
date
