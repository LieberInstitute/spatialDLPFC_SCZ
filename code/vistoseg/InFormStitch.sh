#!/bin/bash
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 5
#$ -o logs/InformStitch_Maddy.txt
#$ -e logs/InformStitch_Maddy.txt
#$ -m e
#$ -M madhavitippani28@gmail.com


echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
# echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id:  V12F14-053_NTC"
echo "****"

## load MATLAB
module load matlab/R2019a

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib/'


## Read parameters
filename='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/round1/20220814_VIF_PNN_S1_NTC/20220814_VIF_PNN_S1_NTC_Scan1_*_component_data.tif'
fname='Channelorder_test'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$filename',O,7,'$fname')"

echo "**** Job ends ****"
date