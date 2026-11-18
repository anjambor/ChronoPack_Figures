#!/usr/bin/env bash
# Steps 01-03 (expression input, BEELINE gene ordering, TFs + N datasets) for one context. Submit from Run_BEELINE/:
#   BEELINE_CONTEXT=THP1 sbatch scripts/slurm/prep.sbatch     (default K562)
#   BEELINE_GENE_FILTER=chronopack ...                        (default beeline; see scripts/common.py)
#SBATCH --job-name=beeline_prep
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Run_BEELINE directory}"
[[ -f run_beeline.sh ]] || { echo "Submit from the Run_BEELINE directory" >&2; exit 2; }
CP_PYTHON="${CP_PYTHON:-/new-stg/home/alex/miniconda3/envs/chronopack/bin/python}"
RSCRIPT="${RSCRIPT:-$PWD/.conda_env/bin/Rscript}"
export BEELINE_CONTEXT="${BEELINE_CONTEXT:-K562}"
export BEELINE_GENE_FILTER="${BEELINE_GENE_FILTER:-beeline}"
echo "context $BEELINE_CONTEXT, gene filter $BEELINE_GENE_FILTER"
PREP_STEPS="${PREP_STEPS:-01 02 03}"  # e.g. PREP_STEPS=03 to rebuild only the TFs + N inputs from existing prep files
echo "steps $PREP_STEPS, BEELINE_PVAL ${BEELINE_PVAL:-0}"
[[ " $PREP_STEPS " == *" 01 "* ]] && { echo "=== 01 expression input $(date '+%F %T')"; "$CP_PYTHON" scripts/01_prepare_expression.py; }
[[ " $PREP_STEPS " == *" 02 "* ]] && { echo "=== 02 gene ordering (gam) $(date '+%F %T')"; "$RSCRIPT" scripts/02_gene_ordering.R; }
[[ " $PREP_STEPS " == *" 03 "* ]] && { echo "=== 03 TFs + N inputs $(date '+%F %T')"; "$CP_PYTHON" scripts/03_generate_inputs.py; }
echo "BEELINE_PREP_OK $BEELINE_CONTEXT $BEELINE_GENE_FILTER p=${BEELINE_PVAL:-0}"
