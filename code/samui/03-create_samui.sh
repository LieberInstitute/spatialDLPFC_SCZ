#!/bin/bash

#SBATCH --mem=20G
#SBATCH --job-name=03-create_samui
#SBATCH -o /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/samui/logs/
#SBATCH -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/samui/logs/
#SBATCH --array=1-63%4

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"s
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

donor=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/samui/samples_list.txt)
echo "Processing sample ${donor}"
date


module load samui/1.0.0-next.45
python 03-create_samui.py $donor

echo "**** Job ends ****"
date

