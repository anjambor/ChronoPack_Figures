#!/usr/bin/env bash
# Step 06: Fig5 04e gene-universe sensitivity with the r9 BEELINE reruns. Submit from Run_BEELINE/:
#   sbatch scripts/slurm/gene_universe.sbatch
#SBATCH --job-name=beeline_gene_universe
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --time=6:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Run_BEELINE directory}"
[[ -f run_beeline.sh ]] || { echo "Submit from the Run_BEELINE directory" >&2; exit 2; }
CP_PYTHON="${CP_PYTHON:-/new-stg/home/alex/miniconda3/envs/chronopack/bin/python}"
"$CP_PYTHON" scripts/06_gene_universe_sensitivity.py
