#!/bin/bash
#$ -l mem_free=100G,h_vmem=100G,h_stack=256M,h_fsize=100G
#$ -o logs/$JOB_NAME.out
#$ -e logs/$JOB_NAME.err
#$ -m e
#$ -M bguo6@jhu.edu
#$ -t 1
##$ -tc 4

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: ${JOB_NAME}"
echo "****"

## load MATLAB
module load matlab/R2019a

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg/lib/'

## Read inputs from refineVNS_list.txt file
fname=$(awk 'BEGIN {FS="\t"} {print $1}' parameters/${JOB_NAME}.tsv | awk "NR==1")
M=$(awk 'BEGIN {FS="\t"} {print $2}' parameters/${JOB_NAME}.tsv | awk "NR==1")

## Run refineVNS function
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), refineVNS('$fname',$M)"

echo "**** Job ends ****"
date