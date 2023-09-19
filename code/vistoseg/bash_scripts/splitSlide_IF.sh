#!/bin/bash
#$ -pe local 10
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o //dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/SS_slide_10_11_12_13_output.txt
#$ -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/SS_slide_10_11_12_13_error.txt
#$ -m be
#$ -M uma.kaipa.94@gmail.com



echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
# echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id:  V13F27-295 V13F27-296 V13M06-340 V13M06-342"
echo "****"

## load MATLAB
module load matlab/R2019a

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/'


## slide 10
fname1='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/11_Real_PNN_Slide10/V13F27-295.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname1')" 
mv /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/11_Real_PNN_Slide10/V13F27-295*_* /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/


## slide 11
fname2='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/12_Real_PNN_Slide11/V13F27-296.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname2')" 
mv /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/12_Real_PNN_Slide11/V13F27-296*_* /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/


## slide 12
fname3='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/13_Real_PNN_Slide12/V13M06-340.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname3')" 
mv /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/13_Real_PNN_Slide12/V13M06-340*_* /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/

## slide 13
fname4='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/14_Real_PNN_Slide13/V13M06-342.mat'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname4')" 
mv /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/14_Real_PNN_Slide13/V13M06-342*_* /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/

