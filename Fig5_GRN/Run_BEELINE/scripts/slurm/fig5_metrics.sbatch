#!/usr/bin/env bash
# Step 05: Fig5 per-TF EPR / AUPRCr of the r9 reruns. Submit from Run_BEELINE/:
#   sbatch scripts/slurm/fig5_metrics.sbatch
#SBATCH --job-name=beeline_fig5_metrics
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Run_BEELINE directory}"
[[ -f run_beeline.sh ]] || { echo "Submit from the Run_BEELINE directory" >&2; exit 2; }
CP_PYTHON="${CP_PYTHON:-/new-stg/home/alex/miniconda3/envs/chronopack/bin/python}"
"$CP_PYTHON" scripts/05_fig5_per_tf_metrics.py
