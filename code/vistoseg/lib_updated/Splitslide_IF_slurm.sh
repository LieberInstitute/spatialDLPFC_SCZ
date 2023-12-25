#!/bin/bash
#SBATCH --mem=50G
#SBATCH --job-name=SplitSlide_slides_15
#SBATCH -o /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/SS_slides_15_output.txt
#SBATCH -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/logs/SS_slides_15_error.txt
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
echo "Sample id:  V13M06-344"
echo "****"

## load MATLAB
module load matlab


## slide 14
## fname1='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/15_Real_PNN_Slide14/V13M06-343.mat'


## slide 15
## fname2='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/16_Real_PNN_Slide15/V13M06-344.mat'

## slide 16
fname3='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/17_Real_PNN_Slide16/20230919_VSPG_PNN_Slide16_Scan2_.mat'


## echo "splitting V13M06-343"
## matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname1')"
## echo "splitting V13M06-344"
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname2')"
echo "splitting V13M06-336"
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname3')"


echo "**** Job ends ****"
date
