#!/usr/bin/env bash
# Step 05: Fig5 composite laid out as ref/FigureComposite.v6.png (expression method set). Submit from Fig5_GRN/:
#   sbatch scripts/slurm/composite.sbatch        (after the pipeline: needs 02 and 02b outputs)
#SBATCH --job-name=fig5v4_composite
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Fig5_GRN directory}"
[[ -f run_fig5_v4.sh ]] || { echo "Submit from the Fig5_GRN directory" >&2; exit 2; }
CP_ENV="${CP_ENV:-/new-stg/home/alex/miniconda3/envs/chronopack}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/output/.mplconfig}"
"$CP_ENV/bin/python" scripts/05_composite_figure.py
