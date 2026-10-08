#!/bin/bash
# Submit gsMap steps 01-12 as a chain of SLURM jobs with dependencies.
# Run from this directory after 00-setup_gsMap_env.sh has finished:
#   bash run_all.sh            # all steps
#   bash run_all.sh 07         # resume from step 07 (earlier outputs exist)
#
# Per-sample array steps (04-08) depend on each other with `aftercorr`, so
# sample N moves to the next step as soon as its own previous task succeeds.

set -eo pipefail
start_step=${1:-01}

submit() {
  # $1 = step number, $2 = dependency string (may be empty), $3 = script
  if [[ "10#$1" -lt "10#${start_step}" ]]; then
    echo ""
    return
  fi
  local dep_flag=()
  [ -n "$2" ] && dep_flag=(--dependency="$2")
  sbatch --parsable "${dep_flag[@]}" "$3"
}

dep() {
  # $1 = type (afterok / aftercorr), remaining = job ids (empty ids skipped)
  local type=$1; shift
  local ids=()
  for id in "$@"; do [ -n "${id}" ] && ids+=("${id}"); done
  if [ ${#ids[@]} -gt 0 ]; then
    echo "${type}:$(IFS=:; echo "${ids[*]}")"
  fi
}

j01=$(submit 01 "" 01-export_spe_to_h5ad.sh)
j02=$(submit 02 "" 02-format_sumstats.sh)
j03=$(submit 03 "$(dep afterok "${j01}")" 03-create_slice_mean.sh)
j04=$(submit 04 "$(dep afterok "${j01}")" 04-find_latent_representations.sh)

# 05 needs this sample's latent rep (aftercorr) and the slice mean (afterok)
dep05=$(dep aftercorr "${j04}")
dep03=$(dep afterok "${j03}")
j05=$(submit 05 "$(IFS=,; tmp=(${dep05} ${dep03}); echo "${tmp[*]}")" 05-latent_to_gene.sh)

j06=$(submit 06 "$(dep aftercorr "${j05}")" 06-generate_ldscore.sh)

dep06=$(dep aftercorr "${j06}")
dep02=$(dep afterok "${j02}")
j07=$(submit 07 "$(IFS=,; tmp=(${dep06} ${dep02}); echo "${tmp[*]}")" 07-spatial_ldsc.sh)

j08=$(submit 08 "$(dep aftercorr "${j07}")" 08-cauchy_combination_per_sample.sh)
j09=$(submit 09 "$(dep afterok "${j07}")" 09-cauchy_combination_across_samples.sh)
j10=$(submit 10 "$(dep afterok "${j08}")" 10-report_rep_samples.sh)
j11=$(submit 11 "$(dep afterok "${j08}" "${j09}")" 11-collect_gsMap_results.sh)
j12=$(submit 12 "$(dep afterok "${j11}")" 12-viz_gsMap_results.sh)

echo "Submitted job ids:"
for s in 01 02 03 04 05 06 07 08 09 10 11 12; do
  v="j${s}"
  echo "  step ${s}: ${!v:-skipped}"
done
