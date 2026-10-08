#!/usr/bin/env bash
## invoke through sbatch with explicit paths; slurm executes a spooled script copy.
set -eo pipefail
if command -v module >/dev/null 2>&1; then module load conda_R/4.5; fi
set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
analysis_dir=$1
project_root=$2
output_dir=$3
workers=${SLURM_CPUS_PER_TASK:-1}
Rscript --vanilla "$analysis_dir/tests/test_core.R"
Rscript --vanilla "$analysis_dir/run_mediation.R" --stage all --project-root "$project_root" --outdir "$output_dir" --workers "$workers"
Rscript --vanilla "$analysis_dir/ptn_by_spd.R" "$project_root" "$output_dir"
python3 "$analysis_dir/tests/verify_outputs.py" --outdir "$output_dir"
