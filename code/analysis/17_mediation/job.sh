#!/usr/bin/env bash
## invoke through sbatch with explicit paths; slurm executes a spooled script copy.
set -eo pipefail
module load conda_R/4.5
set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
analysis_dir=$1
stage=$2
project_root=$3
output_dir=$4
workers=${SLURM_CPUS_PER_TASK:-1}
if [[ "$stage" == "screening" ]]; then
  Rscript --vanilla "$analysis_dir/tests/test_core.R"
  bash "$analysis_dir/tests/test_hit_paths.sh" "$analysis_dir"
  for analysis_stage in audit historical screen; do
    Rscript --vanilla "$analysis_dir/run_mediation.R" --stage "$analysis_stage" --project-root "$project_root" --outdir "$output_dir" --workers "$workers"
  done
else
  Rscript --vanilla "$analysis_dir/run_mediation.R" --stage "$stage" --project-root "$project_root" --outdir "$output_dir" --workers "$workers"
fi

if [[ "$stage" == "report" ]]; then
  python3 "$analysis_dir/tests/verify_outputs.py" --outdir "$output_dir" --require-robustness
fi
