#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=20G
#SBATCH --job-name=create_samui
#SBATCH -o logs/create_samui.txt 
#SBATCH -e logs/create_samui.txt 

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"s
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load samui/1.0.0-next.24
python create_samui.py 

echo "**** Job ends ****"
date
