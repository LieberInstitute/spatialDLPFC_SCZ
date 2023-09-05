#!/bin/bash
#SBATCH --mem=1G
##SBATCH --mem=80G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bguo6@jhu.edu   # TODO: edit this mail address
#SBATCH --output=logs/%x.txt   # file to collect standard output
#SBATCH --error=logs/%x.txt    # file to collect standard output


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_CLUSTER_NAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load SpaceRanger
# TODO: should we use paceranger 2.1 or 2.0
module load spaceranger/2.0.0

## List current modules for reproducibility
module list

## Read parameters
# TODO: change $JOB_NAME to slurm version
# TODO: change {SGE_TASK_ID} to just 1

SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' parameters/${SLURM_JOB_NAME}.tsv)
SLIDE=$(awk 'BEGIN {FS="\t"} {print $2}' parameters/${SLURM_JOB_NAME}.tsv)
CAPTUREAREA=$(awk 'BEGIN {FS="\t"} {print $3}' parameters/${SLURM_JOB_NAME}.tsv)
IMAGEPATH=$(awk 'BEGIN {FS="\t"} {print $4}' parameters/${SLURM_JOB_NAME}.tsv)
LOUPEPATH=$(awk 'BEGIN {FS="\t"} {print $5}' parameters/${SLURM_JOB_NAME}.tsv)
FASTQPATH=$(awk 'BEGIN {FS="\t"} {print $6}' parameters/${SLURM_JOB_NAME}.tsv)

echo "Processing sample ${SAMPLE} from slide ${SLIDE} and capture area ${CAPTUREAREA} with image ${IMAGEPATH} and aligned with ${LOUPEPATH} with FASTQs: ${FASTQPATH}"
date

## For keeping track of dates of the input files
ls -lh ${IMAGEPATH}
ls -lh ${LOUPEPATH}

## Hank from 10x Genomics recommended setting this environment
# export NUMBA_NUM_THREADS=1

## Run SpaceRanger
# TODO: Compare this to Ryan's code
# spaceranger count \
#     --id=${SAMPLE} \
#     --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
#     --fastqs=${FASTQPATH} \
#     --darkimage=${IMAGEPATH} \
#     --slide=${SLIDE} \
#     --area=${CAPTUREAREA} \
#     --loupe-alignment=${LOUPEPATH} \
#     --jobmode=local \
#     --localcores=9 \
#     --localmem=81

## Update file permission
# TODO: uncomment at the end
# sh /dcs04/lieber/lcolladotor/_jhpce_org_LIBD001/update_permissions_spatialteam.sh ${SAMPLE}


## Move output
echo "Moving results to new location"
date
# TODO: check if the path is actually correct becasue moving the folder position
mkdir -p /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/
mv ${SAMPLE} /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/


echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
