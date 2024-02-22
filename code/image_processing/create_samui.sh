#!/bin/bash

#$ -cwd
#$ -N "create_samui"
#$ -o code/image_processing/logs/create_samui.log
#$ -e code/image_processing/logs/create_samui.log
#$ -l mf=80G,h_vmem=80G,h_fsize=50G

#SBATCH -p shared
#SBATCH --mem=80G
#SBATCH --job-name=create_samui
#SBATCH -o code/image_processing/logs/create_samui.log
#SBATCH -e code/image_processing/logs/create_samui.log

if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
    module_name=samui
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
    module_name=loopy
fi


echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

module load ${module_name}/1.0.0-next.24
python create_samui.py 

echo "**** Job ends ****"
date
} > $log_path 2>&1