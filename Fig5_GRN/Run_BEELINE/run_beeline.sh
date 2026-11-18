#!/usr/bin/env bash
# Run_BEELINE: BEELINE (GENIE3, GRNBOOST2, PPCOR, SINCERITIES) on r9 single cells, TFs + 1000 and TFs + 2000 gene
# sets, for each context: K562 (r9 K562 + TNF-alpha) and THP1 (THP1 cells of the r9 THP1 + K562 + LPS experiment).
# Submits Slurm jobs; run from this directory on a login/compute node:
#
#   ./run_beeline.sh                          prep (steps 01-03) per context, then its 8 algorithm jobs after it
#   CONTEXTS=THP1 ./run_beeline.sh            one context
#   SKIP_PREP=1 ./run_beeline.sh              algorithm jobs only (inputs already built)
#   ALGOS="GRNBOOST2 PPCOR" DATASETS="K562-r9-sc-TFs1000" SKIP_PREP=1 ./run_beeline.sh   a subset
#   BEELINE_GENE_FILTER=chronopack ./run_beeline.sh   ChronoPack's expression filter instead of BEELINE's 10% filter
#                                                     (prep in output/prep/<ctx>-cpfilt/, datasets <ctx>-r9-sc-cpfilt-TFs<N>)
#   BEELINE_GENE_FILTER=chronopack BEELINE_PVAL=0.05 PREP_STEPS=03 ./run_beeline.sh
#       the paper's selection with an uncorrected P < 0.05 cut-off (datasets <ctx>-r9-sc-cpfilt-p05-TFs<N>); PREP_STEPS=03
#       reuses the existing prep files (steps 01-02 do not depend on the cut-off)
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output/logs

CONTEXTS="${CONTEXTS:-K562 THP1}"
export BEELINE_GENE_FILTER="${BEELINE_GENE_FILTER:-beeline}"
case "$BEELINE_GENE_FILTER" in beeline) tag="" ;; chronopack) tag="-cpfilt" ;;
  *) echo "BEELINE_GENE_FILTER must be beeline or chronopack" >&2; exit 2 ;; esac
export BEELINE_PVAL="${BEELINE_PVAL:-0}"
tag="$tag$(python3 -c "p=float('$BEELINE_PVAL'); print('' if p == 0 else '-p' + f'{p:g}'.split('.')[-1])")"
ALGOS="${ALGOS:-GENIE3 GRNBOOST2 PPCOR SINCERITIES}"
declare -A RES=(  # sbatch resources per algorithm
  [GENIE3]="-c 40 --mem=90G -t 4-00:00:00 --exclude=bigcompute"  # THP1-TFs2000 hung on bigcompute (2026-10-20)
  [GRNBOOST2]="-c 24 --mem=64G -t 2-00:00:00"
  [PPCOR]="-c 4 --mem=64G -t 1-00:00:00"
  [SINCERITIES]="-c 4 --mem=64G -t 2-00:00:00"
)

for ctx in $CONTEXTS; do
  dep=()
  if [[ "${SKIP_PREP:-0}" != "1" ]]; then
    prep=$(BEELINE_CONTEXT=$ctx sbatch --parsable -J "beeline_prep_$ctx$tag" --export=ALL scripts/slurm/prep.sbatch)
    echo "$ctx prep: job $prep"
    dep=(--dependency=afterok:$prep)
  fi
  datasets="${DATASETS:-$ctx-r9-sc$tag-TFs1000 $ctx-r9-sc$tag-TFs2000}"
  for ds in $datasets; do
    [[ "${ds%%-*}" == "$ctx" ]] || continue  # DATASETS may list several contexts
    for algo in $ALGOS; do
      # shellcheck disable=SC2086
      id=$(DATASET=$ds ALGO=$algo sbatch --parsable "${dep[@]}" -J "bl_${algo}_${ctx}${tag}_${ds##*-}" ${RES[$algo]} \
           --export=ALL scripts/slurm/beeline_algo.sbatch)
      echo "$ds $algo: job $id"
    done
  done
done
