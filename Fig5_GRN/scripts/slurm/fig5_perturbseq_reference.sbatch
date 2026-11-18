#!/usr/bin/env bash
# Submit from Fig5_GRN/:  sbatch scripts/slurm/fig5_perturbseq_reference.sbatch
# Step 02e only (needs the with_prior full and shared-universe tables): figures/perturbseq_reference_v2
#SBATCH --job-name=fig5v4_pseq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Fig5_GRN directory}"
[[ -f run_fig5_v4.sh ]] || { echo "Submit from the Fig5_GRN directory" >&2; exit 2; }
CP_ENV="${CP_ENV:-/new-stg/home/alex/miniconda3/envs/chronopack}"
export PATH="$CP_ENV/bin:$PATH"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/output/.mplconfig}"
export FIG5_METHOD_SET=with_prior
unset FIG5_UNIVERSE
python scripts/02e_perturbseq_reference.py ${PSEQ_ARGS:-}
