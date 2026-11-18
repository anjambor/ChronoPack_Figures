#!/usr/bin/env bash
# Submit from Fig5_GRN/:  sbatch scripts/slurm/fig5_pipeline.sbatch
#SBATCH --job-name=fig5v4_grn
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Fig5_GRN directory}"
[[ -f run_fig5_v4.sh ]] || { echo "Submit from the Fig5_GRN directory" >&2; exit 2; }
./run_fig5_v4.sh
