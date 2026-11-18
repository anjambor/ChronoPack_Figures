#!/usr/bin/env bash
# Submit from Fig3_DynamicModes/:  sbatch scripts/slurm/fig3_v2_pipeline.sbatch
# (fits missing models itself; for a fresh start run the two fit.sbatch jobs first)
#SBATCH --job-name=fig3v2_pipeline
#SBATCH --mem=90G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Fig3_DynamicModes directory}"
[[ -f run_fig3_v2.sh ]] || { echo "Submit from the Fig3_DynamicModes directory" >&2; exit 2; }
./run_fig3_v2.sh
