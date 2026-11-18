#!/usr/bin/env bash
# Submit from FigS_C2C_ligand_specificity/:  sbatch scripts/slurm/ligand_specificity.sbatch
#SBATCH --job-name=c2c_ligand_specificity
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the FigS_C2C_ligand_specificity directory}"
[[ -f run_ligand_specificity.sh ]] || { echo "Submit from the FigS_C2C_ligand_specificity directory" >&2; exit 2; }
./run_ligand_specificity.sh
