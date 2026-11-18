#!/usr/bin/env bash
# Step 11: P < 0.05 round figures from step 10's tables. Submit from Run_BEELINE/: sbatch scripts/slurm/p05_figures.sbatch
#SBATCH --job-name=beeline_p05_figures
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Run_BEELINE directory}"
[[ -f run_beeline.sh ]] || { echo "Submit from the Run_BEELINE directory" >&2; exit 2; }
CP_PYTHON="${CP_PYTHON:-/new-stg/home/alex/miniconda3/envs/chronopack/bin/python}"
"$CP_PYTHON" scripts/11_plot_p05_figures.py
