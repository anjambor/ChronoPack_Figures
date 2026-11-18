#!/usr/bin/env bash
# Submit from Fig6_CrossContext/:  sbatch scripts/slurm/fig6_pipeline.sbatch
#SBATCH --job-name=fig6v2_crossctx
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Fig6_CrossContext directory}"
[[ -f run_fig6.sh ]] || { echo "Submit from the Fig6_CrossContext directory" >&2; exit 2; }
./run_fig6.sh
