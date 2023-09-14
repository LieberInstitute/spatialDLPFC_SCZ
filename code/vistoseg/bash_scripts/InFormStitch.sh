#!/bin/bash
#$ -pe local 5
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o //dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/IS_V13M06-281_output.txt
#$ -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/IS_V13M06-281_error.txt
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
echo "Sample id:  V13M06-281, V13M06-282, V13F27-293, V13F27-294"
echo "****"

## load MATLAB
module load matlab/R2019a

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/'


## Read parameters
path1='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/7_Real_PNN_Slide6/20230717_VSPG_PNN_Slide6_Scan1_*_component_data.tif'
fname1='V13M06-281'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path1',O,6,'$fname1')"

path2='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/8_Real_PNN_Slide7/20230718_VSPG_PNN_Slide7_Scan1_*_component_data.tif'
fname2='V13M06-282'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path2',O,6,'$fname2')"

path3='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/9_Real_PNN_Slide8/20230724_VSPG_PNN_Slide8_Scan1_*_component_data.tif'
fname3='V13F27-293'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path3',O,6,'$fname3')"

path4='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/10_Real_PNN_Slide9/20230725_VSPG_PNN_Slide9_Scan1_*_component_data.tif'
fname4='V13F27-294'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path4',O,6,'$fname4')"


echo "**** Job ends ****"
date