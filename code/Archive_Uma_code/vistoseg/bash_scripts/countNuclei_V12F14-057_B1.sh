#!/bin/bash
#$ -l mem_free=50G,h_vmem=50G,h_fsize=10G
#$ -o /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/CN_V12F14-057_B1_output.txt
#$ -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/CN_V12F14-057_B1_error.txt
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
echo "Sample id: /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/segmented_channels_stitched/Test2/V12F14-057_B1.mat"
echo "****"

## List current modules for reproducibility
module matlab/R2019a
mask='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/segmented_channels_stitched/Test2/V12F14-057_B1.mat'
jsonname='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V12F14-057_B1/outs/spatial/scalefactors_json.json'
posname='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V12F14-057_B1/outs/spatial/tissue_positions.csv'

toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$mask','$jsonname','$posname')" 

echo "**** Job ends ****"
date
