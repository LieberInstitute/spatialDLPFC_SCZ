#!/usr/bin/env bash
## screen first; run independent robustness jobs together; report after both finish.
set -euo pipefail
analysis_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_root=${1:-$(cd "$analysis_dir/../../.." && pwd)}
output_dir=${2:-"$project_root/processed-data/18_mediation"}
mkdir -p "$output_dir/logs"
screen_id=$(sbatch --parsable --job-name=scz-mediation --cpus-per-task=4 --mem=16G --time=04:00:00 \
  --output="$output_dir/logs/screening-%j.log" "$analysis_dir/job.sh" "$analysis_dir" screening "$project_root" "$output_dir")
screen_id=${screen_id%%;*}
sensitivity_id=$(sbatch --parsable --dependency="afterok:$screen_id" --kill-on-invalid-dep=yes \
  --job-name=scz-med-sensitivity --cpus-per-task=4 --mem=16G --time=04:00:00 \
  --output="$output_dir/logs/sensitivity-%j.log" "$analysis_dir/job.sh" "$analysis_dir" sensitivity "$project_root" "$output_dir")
sensitivity_id=${sensitivity_id%%;*}
overlap_id=$(sbatch --parsable --dependency="afterok:$screen_id" --kill-on-invalid-dep=yes \
  --job-name=scz-med-overlap --cpus-per-task=1 --mem=64G --time=04:00:00 \
  --output="$output_dir/logs/overlap-%j.log" "$analysis_dir/job.sh" "$analysis_dir" overlap "$project_root" "$output_dir")
overlap_id=${overlap_id%%;*}
report_id=$(sbatch --parsable --dependency="afterok:$sensitivity_id:$overlap_id" --kill-on-invalid-dep=yes \
  --job-name=scz-med-report --cpus-per-task=1 --mem=8G --time=01:00:00 \
  --output="$output_dir/logs/report-%j.log" "$analysis_dir/job.sh" "$analysis_dir" report "$project_root" "$output_dir")
report_id=${report_id%%;*}
printf 'screening\t%s\nsensitivity\t%s\noverlap\t%s\nreport\t%s\n' \
  "$screen_id" "$sensitivity_id" "$overlap_id" "$report_id" | tee "$output_dir/logs/submitted_jobs.tsv"
