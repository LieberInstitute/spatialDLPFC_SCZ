#!/usr/bin/env bash
## one job: screen, shared-spot refits (loads the 2.6 GB raw spot object), report, PTN by SpD.
set -euo pipefail
analysis_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_root=${1:-$(cd "$analysis_dir/../../.." && pwd)}
output_dir=${2:-"$project_root/processed-data/18_mediation"}
mkdir -p "$output_dir/logs"
sbatch --parsable --job-name=scz-mediation --cpus-per-task=5 --mem=64G --time=02:00:00 \
  --output="$output_dir/logs/mediation-%j.log" "$analysis_dir/job.sh" "$analysis_dir" "$project_root" "$output_dir"
