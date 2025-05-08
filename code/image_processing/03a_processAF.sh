#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=20G
#SBATCH --job-name=AFseg
#SBATCH -o logs/AFseg.txt 
#SBATCH -e logs/AFseg.txt 

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"s
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load MATLAB
module load matlab/R2023a

## Run WFAseg version 2 function
matlab -nodesktop -nosplash -r "run('03a_processAF.m'); exit" 

echo "**** Job ends ****"
date

