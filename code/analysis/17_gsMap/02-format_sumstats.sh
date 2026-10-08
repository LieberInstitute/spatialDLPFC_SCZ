#!/bin/bash
#SBATCH --job-name=02-format_sumstats
#SBATCH --mem=30G
#SBATCH --time=2:00:00
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bguo6@jhu.edu
#SBATCH --output=logs/%x.txt
#SBATCH --error=logs/%x.txt    # file to collect standard output

set -eo pipefail
mkdir -p logs
source gsMap_config.sh
log_job_info
activate_gsmap

mkdir -p "${SUMSTATS_DIR}"

# SCZ (PGC3, European) ----
## Clean raw PGC3 file: drop `##` header, compute N_eff ----
python 02-prep_PGC3_sumstats.py \
  "${PGC3_RAW}" \
  "${SUMSTATS_DIR}/SCZ_PGC3_clean.tsv.gz"

## Convert to gsMap/LDSC format (SNP, A1, A2, Z, N); applies gsMap's default
## INFO >= 0.9 and MAF >= 0.01 filters ----
gsmap format_sumstats \
  --sumstats "${SUMSTATS_DIR}/SCZ_PGC3_clean.tsv.gz" \
  --out "${SUMSTATS_DIR}/SCZ_PGC3" \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --beta BETA \
  --se SE \
  --p P \
  --n N \
  --info INFO \
  --frq FRQ

# error prevention
test -s "${SUMSTATS_DIR}/SCZ_PGC3.sumstats.gz"
zcat "${SUMSTATS_DIR}/SCZ_PGC3.sumstats.gz" | head -n 5

# Additional traits ----
# Format any extra GWAS listed in traits.tsv here, e.g. a height negative
# control:
# gsmap format_sumstats --sumstats <raw_height_file> --out "${SUMSTATS_DIR}/Height"

log_job_end
