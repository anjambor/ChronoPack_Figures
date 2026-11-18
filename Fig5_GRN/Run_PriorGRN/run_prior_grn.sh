#!/usr/bin/env bash
# Run_PriorGRN driver (run from Run_PriorGRN/ on the login node; it only submits Slurm jobs):
#   ./run_prior_grn.sh setup                 one-time: containers, cisTarget files, .conda_env (setup.sbatch)
#   ./run_prior_grn.sh all [METHODS...]      step 01, then each method on each dataset (default: all three)
#   ./run_prior_grn.sh methods [METHODS...]  methods only (inputs already prepared)
#   AFTER=<jobid> ./run_prior_grn.sh all      start step 01 only after job <jobid> succeeds (e.g. the setup job)
# Datasets: {K562,THP1}-r9-sc-cpfilt-p05-TFs{1000,2000} (scripts/common.py).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output/logs
mode=${1:-all}; shift || true
methods=("${@:-PYSCENIC CELLORACLE INFERELATOR}"); methods=(${methods[*]})
datasets=(K562-r9-sc-cpfilt-p05-TFs1000 K562-r9-sc-cpfilt-p05-TFs2000 THP1-r9-sc-cpfilt-p05-TFs1000 THP1-r9-sc-cpfilt-p05-TFs2000)

if [[ $mode == setup ]]; then sbatch scripts/slurm/setup.sbatch; exit; fi
dep=()
if [[ $mode == all ]]; then
  p=$(sbatch --parsable ${AFTER:+--dependency=afterok:$AFTER} scripts/slurm/prep.sbatch); echo "prep $p"
  dep=(--dependency=afterok:$p)
fi
for m in "${methods[@]}"; do
  case $m in
    PYSCENIC)    res=(-c 32 --mem=96G -t 3-00:00:00) ;;
    CELLORACLE)  res=(-c 16 --mem=96G -t 2-00:00:00) ;;
    INFERELATOR) res=(-c 32 --mem=96G -t 3-00:00:00) ;;
    *) echo "unknown method $m" >&2; exit 2 ;;
  esac
  for ds in "${datasets[@]}"; do
    j=$(METHOD=$m DATASET=$ds sbatch --parsable "${dep[@]}" "${res[@]}" -J "pg_${m}_${ds%%-*}_${ds##*-}" \
        --export=ALL,METHOD=$m,DATASET=$ds scripts/slurm/method.sbatch)
    echo "$m $ds $j"
  done
done
