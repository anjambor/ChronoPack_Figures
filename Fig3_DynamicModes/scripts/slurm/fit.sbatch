#!/usr/bin/env bash
# One ChronoPack downsampling fit. Submit from Fig3_DynamicModes/:
#   RUN=K562_trunc sbatch scripts/slurm/fit.sbatch     # reference, ~1-2 h
#   RUN=K562_ds    sbatch scripts/slurm/fit.sbatch     # downsampled, minutes
#SBATCH --job-name=fig3v2_fit
#SBATCH --mem=90G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Fig3_DynamicModes directory}"
[[ -f run_fig3_v2.sh || -d scripts ]] || { echo "Submit from the Fig3_DynamicModes directory" >&2; exit 2; }
CP_ENV="${CP_ENV:-/new-stg/home/alex/miniconda3/envs/chronopack}"
export PATH="$CP_ENV/bin:$PATH"  # R for rpy2
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/output/.mplconfig}"
"$CP_ENV/bin/python" scripts/01_fit_downsampling.py --run "${RUN:?set RUN=K562_trunc or K562_ds}"
