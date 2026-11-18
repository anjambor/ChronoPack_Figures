#!/usr/bin/env bash
# Submit from Fig4_C2C/:  sbatch scripts/slurm/fig4_v2_pipeline.sbatch   (~5 min)
#SBATCH --job-name=fig4v2_pipeline
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Fig4_C2C directory}"
[[ -f run_fig4_v2.sh ]] || { echo "Submit from the Fig4_C2C directory" >&2; exit 2; }
./run_fig4_v2.sh
