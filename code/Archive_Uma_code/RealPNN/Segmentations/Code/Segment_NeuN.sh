#!/bin/bash
#$ -l caracol,mem_free=120G,h_vmem=120G,h_fsize=10G
#$ -o /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/RealPNN/Segmentations/logs/Output_NeuN_seg.txt
#$ -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/RealPNN/Segmentations/logs/Error_NeuN_seg.txt
#$ -m be
#$ -M uma.kaipa.94@gmail.com


 
echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"

### load the python environment
eval "$(conda shell.bash hook)"
conda activate venv_segmentations

### python filename
python NeuN.py
