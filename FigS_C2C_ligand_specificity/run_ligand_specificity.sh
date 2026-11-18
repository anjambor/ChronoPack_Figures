#!/usr/bin/env bash
# FigS_C2C_ligand_specificity: ChronoPack vs best-comparator ligand ranks on Fig. 4's r9 rankings.
# Run on a compute node:  sbatch scripts/slurm/ligand_specificity.sbatch   (~1 min)
set -euo pipefail
cd "$(dirname "$0")"
PY="${CP_PYTHON:-$HOME/miniconda3/envs/chronopack/bin/python}"
mkdir -p output/logs
"$PY" scripts/01_ligand_specificity.py
"$PY" scripts/02_compare_to_capsule.py
echo LIGAND_SPECIFICITY_PIPELINE_OK
