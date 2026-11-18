#!/usr/bin/env bash
# Step 01 (per-dataset inputs for all four datasets). Submit from Run_PriorGRN/: sbatch scripts/slurm/prep.sbatch
#SBATCH --job-name=priorgrn_prep
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Run_PriorGRN directory}"
[[ -f run_prior_grn.sh ]] || { echo "Submit from the Run_PriorGRN directory" >&2; exit 2; }
export PYTHONNOUSERSITE=1
.conda_env/bin/python scripts/01_prepare_inputs.py
