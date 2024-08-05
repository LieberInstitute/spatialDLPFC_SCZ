#!/bin/bash
#SBATCH --mem=30G
#SBATCH --job-name=countNuclei
#SBATCH -o logs/countNuclei_%a.txt 
#SBATCH -e logs/countNuclei_%a.txt
#SBATCH --array=24
#SBATCH --time=48:00:00
#SBATCH --mail-type=FAIL            # Send email on job failure
#SBATCH --mail-user=madhavitippani28@gmail.com  # Your email address



echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"s
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load MATLAB
module load matlab/R2023a

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/image_processing/code/'

filePath=$(ls -1 /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/DAPI_NeuN_WFA_Segs/*_segs.mat | sed -n "${SLURM_ARRAY_TASK_ID}p")

## Check if file path exists
if [ ! -f "$filePath" ]; then
    echo "File not found: $filePath"
    exit 1
fi

fileName=$(basename "$filePath" _segs.mat)
echo "Processing sample ${fileName}"


## Read inputs
jsonname1=/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/${fileName}/outs/spatial/scalefactors_json.json
posname1=/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/${fileName}/outs/spatial/tissue_spot_counts.csv
imgname=/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/${fileName}.mat

## Run refineVNS function
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$filePath','$imgname','$jsonname1','$posname1')" 

echo "**** Job ends ****"
date



