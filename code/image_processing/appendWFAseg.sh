#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=40G
#SBATCH --job-name=WFAseg
#SBATCH -o logs/WFAseg1.txt 
#SBATCH -e logs/WFAseg1.txt 

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
matlab -nodesktop -nosplash -r "run('appendWFAseg.m'); exit" 

echo "**** Job ends ****"
date



