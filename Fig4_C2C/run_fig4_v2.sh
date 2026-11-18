#!/usr/bin/env bash
# Fig4 v2: Fig4_C2C with the trajectory panels drawn from ChronoPack NB-GAM fits (r9) instead of
# per-time-point means of the log1p layer. Rankings, panels A-B and the manifest are v1's.
#
#   01 ChronoPack fits/curves (THP1 ligands: quick fit; K562 NFKBIA: stored r9 fit); skipped if
#      output/chronopack/fitted_curves.csv exists, unless REFIT=1
#   02 figure   03 comparison with v1
#
# Run on a compute node: sbatch scripts/slurm/fig4_v2_pipeline.sbatch (from this directory).
set -euo pipefail

cd "$(dirname "$0")"

CP_ENV="${CP_ENV:-/new-stg/home/alex/miniconda3/envs/chronopack}"
CP_PYTHON="$CP_ENV/bin/python"
export PATH="$CP_ENV/bin:$PATH"  # R for rpy2
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/output/.mplconfig}"

step() { echo; echo "=== $* ($(date '+%F %T')) ==="; }

if [[ "${REFIT:-0}" == "1" || ! -f output/chronopack/fitted_curves.csv ]]; then
  step "01 ChronoPack fitted curves"
  "$CP_PYTHON" scripts/01_fit_curves.py
else
  echo "Using existing output/chronopack/fitted_curves.csv (REFIT=1 to refit)"
fi

step "02 figure"
"$CP_PYTHON" scripts/02_plot_fig4_v2.py

step "03 comparison with v1"
"$CP_PYTHON" scripts/03_compare_to_v1.py

echo "FIG4_V2_PIPELINE_OK"
