#!/bin/bash
# #$ -pe local 5
#$ -l mem_free=120G,h_vmem=120G,h_fsize=100G
#$ -o /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/logs/Informstitch_output.txt
#$ -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/logs/Informstitch_error.txt
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

module load matlab/R2019a

toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib/'

path1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/round1/20220814_VIF_PNN_S1_NTC/20220814_VIF_PNN_S1_NTC_Scan1_*_component_data.tif'
fname = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/round1/20220814_VIF_PNN_S1_NTC/20220814_VIF_PNN_S1_NTC_Scan1_[6502,29357]_component_data.tif'
O1 = extractMD(fname)
filename1  = 'V12F14-053_test'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), InFormStitch('$path1',O1,7,'$filename1')"

echo "**** Job ends ****"
date

