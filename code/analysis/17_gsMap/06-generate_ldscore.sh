#!/bin/bash
#SBATCH --job-name=06-generate_ldscore
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --array=1-63%15
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x_%a.txt
#SBATCH --error=logs/%x_%a.txt    # file to collect standard output

# Per sample: map GSS from genes to SNPs (TSS +/- GENE_WINDOW_SIZE, using
# gsMap's "TSS only" strategy) and compute spot-level stratified LD scores
# against the 1000G EUR reference, for chromosomes 1-22.
# This is the most compute-heavy step; tune --mem/--time after a test sample.
# Requires: 05 for this sample.
# Output: {WORKDIR}/{sample}/generate_ldscore/

set -eo pipefail
mkdir -p logs
source gsMap_config.sh
log_job_info
activate_gsmap

sample_id=$(get_sample_id "${SLURM_ARRAY_TASK_ID}")
echo "Sample: ${sample_id}"

for chrom in $(seq 1 22); do
  echo "---- chromosome ${chrom} ----"
  date
  gsmap run_generate_ldscore \
    --workdir "${WORKDIR}" \
    --sample_name "${sample_id}" \
    --chrom "${chrom}" \
    --bfile_root "${BFILE_ROOT}" \
    --keep_snp_root "${KEEP_SNP_ROOT}" \
    --gtf_annotation_file "${GTF_FILE}" \
    --gene_window_size "${GENE_WINDOW_SIZE}"
done

log_job_end
