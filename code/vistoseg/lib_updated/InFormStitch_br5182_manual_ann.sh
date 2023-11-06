#!/bin/bash
#SBATCH --mem=20G
#SBATCH --job-name=InFormStitch_br5182_manual_annoataion_single_section_stitch
#SBATCH -o /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/IS_br5182.log
#SBATCH -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/IS_br5182.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=uma.kaipa.94@gmail.com
#SBATCH -c 5

if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
fi

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"
echo "****"
echo "Sample id:  br5182 br2039"
echo "****"

## load MATLAB
module load matlab

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/'


## Read parameters

## Br 5182
path1='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining/20220712_VIF_MockPNN_Strong_NTC_Scan1_*_component_data.tif'
fname1='V12F14-053'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path1',O,6,'$fname1')"


## Br 2039
path2='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/20220712_VIF_MockPNN_Strong_SCZ_C1_Br2039_MLtraining/20220712_VIF_MockPNN_Strong_SCZ_Scan1_*_component_data.tif'
fname2='V12F14-057'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path2',O,6,'$fname2')"

