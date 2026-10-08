#!/bin/bash
# Shared configuration for the 17_gsMap pipeline.
# Sourced by every step's .sh script (`source gsMap_config.sh`); not run directly.
# Edit paths/parameters here only, so all steps stay consistent.

# Paths ----
REPO_ROOT=$(git rev-parse --show-toplevel)
CODE_DIR="${REPO_ROOT}/code/analysis/17_gsMap"

# Large intermediate files (h5ad, LD scores, resources) live here and are
# git-ignored. Small, summarized R outputs go to processed-data/rds/17_gsMap.
GSMAP_DIR="${REPO_ROOT}/processed-data/17_gsMap"
H5AD_DIR="${GSMAP_DIR}/h5ad"             # per-sample h5ad (01)
SUMSTATS_DIR="${GSMAP_DIR}/sumstats"     # gsMap-formatted GWAS (02)
SLICE_MEAN_FILE="${GSMAP_DIR}/slice_mean/spe_slice_mean.parquet" # (03)
WORKDIR="${GSMAP_DIR}/workdir"           # gsMap --workdir (04-08, 10)
CAUCHY_ACROSS_DIR="${GSMAP_DIR}/cauchy_across_samples" # (09)

SAMPLE_LIST="${H5AD_DIR}/sample_list.txt"         # one sample_id per line
SAMPLE_LIST_NTC="${H5AD_DIR}/sample_list_ntc.txt"
SAMPLE_LIST_SCZ="${H5AD_DIR}/sample_list_scz.txt"
TRAITS_FILE="${CODE_DIR}/traits.tsv"              # trait_name <TAB> sumstats file

# Raw PGC3 SCZ GWAS (hg19), same file used in 15_eqtl_coloc
PGC3_RAW="${REPO_ROOT}/processed-data/ref/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"

# gsMap reference resources (downloaded in 00-setup_gsMap_env.sh) ----
RESOURCE_DIR="${GSMAP_DIR}/gsMap_resource"
BFILE_ROOT="${RESOURCE_DIR}/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC"
KEEP_SNP_ROOT="${RESOURCE_DIR}/LDSC_resource/hapmap3_snps/hm"
GTF_FILE="${RESOURCE_DIR}/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf"
W_FILE="${RESOURCE_DIR}/LDSC_resource/weights_hm3_no_hla/weights."

# gsMap parameters ----
GSMAP_VERSION="1.73.8"
GSMAP_ENV="gsMap"            # conda env name (or full prefix path)
ANNOTATION="spd_label"       # PRECAST-07 spatial domains, obs column in h5ad
DATA_LAYER="count"           # raw UMI counts layer in h5ad
NUM_NEIGHBOUR=51             # latent-space neighbours (gsMap default)
NUM_NEIGHBOUR_SPATIAL=201    # spatial neighbours (gsMap default)
GENE_WINDOW_SIZE=50000       # TSS +/- 50kb SNP-to-gene window (gsMap default)

# Representative samples for gsMap HTML reports (Br8667 NTC, Br5973 SCZ)
REP_SAMPLES=("V13M06-342_D1" "V13M06-343_D1")


# Helper functions ----
log_job_info() {
  echo "**** Job starts ****"
  date
  echo "**** JHPCE info ****"
  echo "User: ${USER}"
  echo "Job id: ${SLURM_JOBID}"
  echo "Job name: ${SLURM_JOB_NAME}"
  echo "Hostname: ${SLURM_CLUSTER_NAME}"
  echo "Task id: ${SLURM_ARRAY_TASK_ID}"
}

activate_gsmap() {
  module load conda
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${GSMAP_ENV}"
  echo "gsMap version: $(python -c 'import gsMap; print(gsMap.__version__)')"
}

# Sample id on line N of the sample list (N = SLURM_ARRAY_TASK_ID)
get_sample_id() {
  local sample_id
  sample_id=$(sed -n "${1}p" "${SAMPLE_LIST}")
  if [ -z "${sample_id}" ]; then
    echo "No sample on line ${1} of ${SAMPLE_LIST}" >&2
    exit 1
  fi
  echo "${sample_id}"
}

# Trait names from traits.tsv (skip comments / header / blank lines)
get_trait_names() {
  grep -v -e '^#' -e '^trait_name' -e '^[[:space:]]*$' "${TRAITS_FILE}" | cut -f1
}

get_sumstats_file() {
  local file
  file=$(grep -v '^#' "${TRAITS_FILE}" | awk -F'\t' -v t="$1" '$1 == t {print $2}')
  echo "${SUMSTATS_DIR}/${file}"
}

log_job_end() {
  echo "**** Job ends ****"
  date
}
