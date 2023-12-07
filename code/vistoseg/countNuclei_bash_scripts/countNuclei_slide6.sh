#!/bin/bash
#SBATCH --mem=50G
#SBATCH --job-name=countNuclei_slide6
#SBATCH -o /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/countNuclei_bash_scripts/logs/CN_slides_6_output.txt
#SBATCH -e /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/countNuclei_bash_scripts/logs/CN_slides_6_error.txt
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
echo "Sample id: V13M06-281"
echo "****"

## load MATLAB
module load matlab

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/'

## tissue position A1
mask1='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/segmented_channels_stitched/slide6/V13M06-281_A1_thresholded.mat'
jsonname1='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V13M06-281_A1/outs/spatial/scalefactors_json.json'
posname1='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V13M06-281_A1/outs/spatial/tissue_positions.csv'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei_without_WFA('$mask1','$jsonname1','$posname1')" 


## tissue position B1
mask2='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/segmented_channels_stitched/slide6/V13M06-281_B1_thresholded.mat'
jsonname2='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V13M06-281_B1/outs/spatial/scalefactors_json.json'
posname2='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V13M06-281_B1/outs/spatial/tissue_positions.csv'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei_without_WFA('$mask2','$jsonname2','$posname2')" 


## tissue position C1
mask3='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/segmented_channels_stitched/slide6/V13M06-281_C1_thresholded.mat'
jsonname3='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V13M06-281_C1/outs/spatial/scalefactors_json.json'
posname3='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V13M06-281_C1/outs/spatial/tissue_positions.csv'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei_without_WFA('$mask3','$jsonname3','$posname3')" 


## tissue position D1
mask4='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/segmented_channels_stitched/slide6/V13M06-281_D1_thresholded.mat'
jsonname4='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V13M06-281_D1/outs/spatial/scalefactors_json.json'
posname4='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V13M06-281_D1/outs/spatial/tissue_positions.csv'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei_without_WFA('$mask4','$jsonname4','$posname4')" 


echo "**** Job ends ****"
date
