#!/bin/bash
#$ -l mem_free=100G,h_vmem=100G,h_fsize=10G
#$ -o /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/logs/countSpots_output_V12F053_A1_centroid.txt
#$ -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/logs/countSpots_error_V12F053_A1_centroid.txt
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
echo "Sample id: /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/all_channels_segemented/Test2/V12F14-053_A1.mat"
echo "****"

## List current modules for reproducibility
module matlab/R2019a
mask='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/all_channels_segemented/Test2/V12F14-053_A1.mat'
jsonname='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V12F14-053_A1/outs/spatial/scalefactors_json.json'
posname='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V12F14-053_A1/outs/spatial/tissue_positions.csv'

toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_centroid/'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$mask','$jsonname','$posname')" 

echo "**** Job ends ****"
date
