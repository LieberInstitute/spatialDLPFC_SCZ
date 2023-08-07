#!/bin/bash
#$ -pe local 5
#$ -l mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -o //dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/IS_V13M06-279_output.txt
#$ -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/IS_V13M06-279_error.txt
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
echo "Sample id:  V13M06-279"
echo "****"

## load MATLAB
module load matlab/R2019a

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated'


## Read parameters
path1='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/5_Real_PNN_Slide4/20230531_VSPG_PNN_Round2_Scan1_[5979,27982]_component_data.tif'

fname1='V13M06-279'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path1',O,6,'$fname1')"

echo "**** Job ends ****"
date