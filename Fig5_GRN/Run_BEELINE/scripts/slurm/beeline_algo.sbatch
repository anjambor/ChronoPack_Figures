#!/usr/bin/env bash
# One BEELINE algorithm on one dataset (BLRunner.py with its single-algorithm config). Submit from Run_BEELINE/:
#   DATASET=K562-r9-sc-TFs1000 ALGO=GRNBOOST2 sbatch scripts/slurm/beeline_algo.sbatch
# run_beeline.sh sets per-algorithm resources with sbatch -c / --mem / -t.
#SBATCH --job-name=beeline_algo
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Run_BEELINE directory}"
[[ -f run_beeline.sh ]] || { echo "Submit from the Run_BEELINE directory" >&2; exit 2; }
: "${DATASET:?}" "${ALGO:?}"
BEELINE_ENV="${BEELINE_ENV:-/new-stg/home/alex/miniconda3/envs/BEELINE}"
export PATH="$BEELINE_ENV/bin:$PATH"
mkdir -p beeline/tmp
export TMPDIR="$PWD/beeline/tmp"  # BLRunner's per-run temp dirs (bound to /tmp in the containers)
unset XDG_RUNTIME_DIR
cd beeline/Beeline
export PYTHONPATH="${PYTHONPATH:-}:$PWD"
echo "BEELINE $ALGO on $DATASET: job $SLURM_JOB_ID on $(hostname), $SLURM_CPUS_PER_TASK CPUs, $(date '+%F %T')"
python BLRunner.py --config "../BEELINE-data/config-files/scRNA-Seq/$DATASET/$ALGO.yaml"
out="../outputs/scRNA-Seq/$DATASET/$ALGO"
ls -la "$out"
[[ -s "$out/rankedEdges.csv" ]] || { echo "missing $out/rankedEdges.csv" >&2; exit 1; }
echo "BEELINE_ALGO_OK $DATASET $ALGO $(date '+%F %T')"
