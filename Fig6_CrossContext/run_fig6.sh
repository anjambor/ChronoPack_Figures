#!/usr/bin/env bash
# Fig6_CrossContext_v2: Fig6_CrossContext with the regression-screen panels drawn from the v4 kinetics tables of
# Fig6_Reproduction/reference/ (superseding Run3), plus the per-family (family-only) CV R^2 figure.
# ORA + covariate + PLS-motif dot plot (v2 percentile covariates, q < 0.05 outlines), panels 6d/6e/6de, the
# full-model CV R^2 figure and the family-only CV R^2 figure. Inputs are the frozen reference tables copied into input/.
#
# Run on a compute node: sbatch scripts/slurm/fig6_pipeline.sbatch (from this directory).
set -euo pipefail

cd "$(dirname "$0")"

PYTHON="${PYTHON:-/new-stg/home/alex/miniconda3/envs/dynamic_mode_motifs/bin/python}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/output/.mplconfig}"

step() { echo; echo "=== $* ($(date '+%F %T')) ==="; }

step "ORA + covariate + PLS motif dot plot (v2, q05)"
"$PYTHON" scripts/plot_fig6_dotplot.py --variant v2 --cutoffs q05

step "regression screen panels d, e, d+e composite (v4 tables)"
"$PYTHON" scripts/plot_fig6de.py --panel d e de

step "full-model CV R^2 per outcome (publication version of refs/figureA_full_model_cv_r2.png)"
"$PYTHON" scripts/plot_fig6_full_model_r2.py

step "family-only CV R^2 per family x outcome (gain over cell type alone; raw)"
"$PYTHON" scripts/plot_fig6_family_only_r2.py

step "checks against the reference manifest and Fig6_CrossContext (v1)"
"$PYTHON" scripts/compare_to_v1.py

echo "FIG6_V2_PIPELINE_OK"
