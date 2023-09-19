#!/bin/bash
#$ -pe local 5
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o //dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/IS_slides_10_11_12_13_output.txt
#$ -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/IS_slides_10_11_12_13_error.txt
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


## Read parameters

## slide 10
path1='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/11_Real_PNN_Slide10/20230829_VSPG_PNN_Slide10_Scan1_*_component_data.tif'
fname1='V13F27-295'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path1',O,6,'$fname1')"

## slide 11
path2='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/12_Real_PNN_Slide11/20230830_VSPG_PNN_Slide11_Scan1_*_component_data.tif'
fname2='V13F27-296'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path2',O,6,'$fname2')"

## slide 12
path3='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/13_Real_PNN_Slide12/20230905_VSPG_PNN_Slide12_Scan1_*_component_data.tif'
fname3='V13M06-340'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path3',O,6,'$fname3')"


## slide 13
path4='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/14_Real_PNN_Slide13/20230906_VSPG_PNN_Slide13_Scan1_*_component_data.tif'
fname4='V13M06-342'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path4',O,6,'$fname4')"


echo "**** Job ends ****"
date