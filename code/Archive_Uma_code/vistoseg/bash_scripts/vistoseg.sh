#!/bin/bash
#$ -l mem_free=100G,h_vmem=100G,h_stack=256M,h_fsize=100G
#$ -o logs/$JOB_NAME.txt
#$ -e logs/$JOB_NAME.txt
#$ -m e
#$ -M bguo6@jhu.edu
##$ -t 1
##$ -tc 4

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
# echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: ${JOB_NAME}"
echo "****"

## load MATLAB
module load matlab/R2019a

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib/'


## Read parameters
mask=$(awk 'BEGIN {FS="\t"} {print $1}' parameters/${JOB_NAME}.tsv | awk "NR==1")
jsonname=$(awk 'BEGIN {FS="\t"} {print $2}' parameters/${JOB_NAME}.tsv | awk "NR==1")
posname=$(awk 'BEGIN {FS="\t"} {print $3}' parameters/${JOB_NAME}.tsv | awk "NR==1")


matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$mask','$jsonname', '$posname')"

echo "**** Job ends ****"
date