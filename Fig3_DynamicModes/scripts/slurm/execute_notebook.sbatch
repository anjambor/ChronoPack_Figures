#!/usr/bin/env bash
# Execute Fig3_v2.ipynb in place on a compute node (chronopack kernel), from Fig3_DynamicModes/:
#   sbatch scripts/slurm/execute_notebook.sbatch                 # show existing outputs
#   FIG_RERUN=1 sbatch scripts/slurm/execute_notebook.sbatch     # regenerate every step
# NB_OUTPUT_DIR=<dir> writes the executed copy there instead of overwriting Fig3_v2.ipynb.
#SBATCH --job-name=fig3v2_notebook
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Fig3_DynamicModes directory}"
[[ -f Fig3_v2.ipynb ]] || { echo "Submit from the Fig3_DynamicModes directory" >&2; exit 2; }
JUPYTER="${JUPYTER:-/new-stg/home/alex/miniconda3/envs/tradeseq_ce/bin/jupyter}"  # nbconvert lives here
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/output/.mplconfig}"
if [[ -n "${NB_OUTPUT_DIR:-}" ]]; then target=(--output-dir "$NB_OUTPUT_DIR"); else target=(--inplace); fi
"$JUPYTER" nbconvert --to notebook --execute "${target[@]}" Fig3_v2.ipynb \
  --ExecutePreprocessor.kernel_name=chronopack --ExecutePreprocessor.timeout=-1
echo "NOTEBOOK_OK FIG_RERUN=${FIG_RERUN:-0}"
