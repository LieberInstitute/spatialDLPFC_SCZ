#!/bin/bash
#SBATCH --job-name=07-spatial_ldsc
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --array=1-63%20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x_%a.txt
#SBATCH --error=logs/%x_%a.txt    # file to collect standard output

# Per sample x trait: spatial stratified LDSC, testing each spot for
# enrichment of GWAS heritability in genes specifically expressed in it.
# Loops over all traits in traits.tsv.
# Requires: 02 (sumstats) and 06 for this sample.
# Output: {WORKDIR}/{sample}/spatial_ldsc/{sample}_{trait}.csv.gz
#         (one row per spot: spot, beta, se, z, p, ...)

set -eo pipefail
mkdir -p logs
source gsMap_config.sh
log_job_info
activate_gsmap

sample_id=$(get_sample_id "${SLURM_ARRAY_TASK_ID}")
echo "Sample: ${sample_id}"

for trait in $(get_trait_names); do
  sumstats_file=$(get_sumstats_file "${trait}")
  echo "---- trait ${trait}: ${sumstats_file} ----"
  date

  # error prevention
  test -s "${sumstats_file}"

  gsmap run_spatial_ldsc \
    --workdir "${WORKDIR}" \
    --sample_name "${sample_id}" \
    --trait_name "${trait}" \
    --sumstats_file "${sumstats_file}" \
    --w_file "${W_FILE}" \
    --num_processes "${SLURM_CPUS_PER_TASK}"
done

log_job_end
