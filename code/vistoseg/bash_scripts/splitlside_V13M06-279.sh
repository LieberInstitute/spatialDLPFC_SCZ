#!/bin/bash
#$ -l mem_free=80G,h_vmem=80G,h_fsize=100G
#$ -o /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/SS_V13M06-279_output.txt
#$ -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/SS_V13M06-279_error.txt
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
echo "Sample id: /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/5_Real_PNN_Slide4/V13M06_279.mat"
echo "****"

## List current modules for reproducibility
module matlab/R2019a
fname='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/5_Real_PNN_Slide4/V13M06_279.mat'

toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname')" 