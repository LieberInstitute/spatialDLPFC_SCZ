#!/bin/bash
#SBATCH --mem=40G
#SBATCH --job-name=InFormStitch_slides_15_16
#SBATCH -o /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/IS_slides_15_16_output.txt
#SBATCH -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/IS_slides_15_16_error.txt
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
echo "Sample id: V13M06-344 V13F27-336"
echo "****"

## load MATLAB
module load matlab

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/'


## Read parameters

## slide 14
## path1='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/15_Real_PNN_Slide14/20230911_VSPG_PNN_Slide14_Scan1_*_component_data.tif'
## fname1='V13M06-343'
## matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path1',O,6,'$fname1')"

## slide 15
path2='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/16_Real_PNN_Slide15/20230913_VSPG_PNN_Slide15_Scan1_*_component_data.tif'
fname2='V13M06-344'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path2',O,6,'$fname2')"

## slide 16
path3='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/17_Real_PNN_Slide16/20230919_VSPG_PNN_Slide16_Scan2_*_component_data.tif'
fname3='V13F27-336'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'Claudin5'; O{3} = 'NeuN'; O{4} = 'WFA'; O{5} = 'AF'; InFormStitch('$path3',O,6,'$fname3')"


echo "**** Job ends ****"
date
