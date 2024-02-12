#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o logs/countNuclei.txt 
#SBATCH -e logs/countNuclei.txt
#SBATCH --array=1

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
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib_updated/'

filePath=$(ls -1 /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/image_histograms/*_WFAseg.mat | sed -n "${SLURM_ARRAY_TASK_ID}p")
fileName=$(basename "$filePath" _WFAseg.mat)
echo "Processing sample ${fileName}"

fname="${SAMPLE}_WFAseg.mat"

## Read inputs
jsonname1=/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/${fileName}/outs/spatial/scalefactors_json.json
posname1=/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/${fileName}/outs/spatial/tissue_positions.csv


## Run refineVNS function
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$filePath','$jsonname1','$posname1')" 

echo "**** Job ends ****"
date



