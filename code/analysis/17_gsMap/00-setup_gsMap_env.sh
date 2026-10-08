#!/bin/bash
# One-time setup: create the gsMap conda env and download gsMap reference
# resources. Run interactively on a JHPCE compute node (not a login node):
#   srun --pty --mem=10G --time=4:00:00 bash
#   bash 00-setup_gsMap_env.sh

set -eo pipefail
source gsMap_config.sh

# Create conda env ----
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"

if conda env list | grep -qE "^${GSMAP_ENV}[[:space:]]"; then
  echo "conda env '${GSMAP_ENV}' already exists, skipping creation"
else
  conda create -y -n "${GSMAP_ENV}" python=3.11
fi

conda activate "${GSMAP_ENV}"
pip install "gsMap==${GSMAP_VERSION}"
gsmap --help | head -n 20

## Record exact environment for reproducibility ----
mkdir -p logs
pip freeze > logs/gsMap_env_pip_freeze.txt


# Download gsMap resources (~ tens of GB: 1000G EUR plink, HapMap3 SNPs,
# LDSC weights, GTF, enhancer annotations) ----
mkdir -p "${GSMAP_DIR}"
if [ -d "${RESOURCE_DIR}" ]; then
  echo "${RESOURCE_DIR} already exists, skipping download"
else
  cd "${GSMAP_DIR}"
  wget -c https://yanglab.westlake.edu.cn/data/gsMap/gsMap_resource.tar.gz
  tar -xvzf gsMap_resource.tar.gz
  rm gsMap_resource.tar.gz
  cd "${CODE_DIR}"
fi

## Error prevention: resource files referenced in gsMap_config.sh exist ----
for f in "${BFILE_ROOT}.1.bed" "${KEEP_SNP_ROOT}.1.snp" "${GTF_FILE}" "${W_FILE}1.l2.ldscore.gz"; do
  if [ ! -e "${f}" ]; then
    echo "Missing expected resource file: ${f}" >&2
    echo "Check the layout of ${RESOURCE_DIR} and update gsMap_config.sh" >&2
    exit 1
  fi
done

echo "gsMap environment and resources are ready"
