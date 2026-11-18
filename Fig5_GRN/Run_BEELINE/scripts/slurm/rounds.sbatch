#!/usr/bin/env bash
# Step 10: Fig5 5b metrics + 04e matched-gene analysis for every BEELINE round, K562 and THP1. Submit from Run_BEELINE/:
#   sbatch scripts/slurm/rounds.sbatch
#SBATCH --job-name=beeline_rounds
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Run_BEELINE directory}"
[[ -f run_beeline.sh ]] || { echo "Submit from the Run_BEELINE directory" >&2; exit 2; }
CP_PYTHON="${CP_PYTHON:-/new-stg/home/alex/miniconda3/envs/chronopack/bin/python}"
"$CP_PYTHON" scripts/10_score_all_rounds.py
