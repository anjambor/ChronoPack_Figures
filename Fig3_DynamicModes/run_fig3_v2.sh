#!/usr/bin/env bash
# Fig3 v2: K562 + TNF-alpha downsampling timing composite rebuilt with ChronoPack fits on r9 data.
#
#   01 fits (reference 0-90 min / 8 knots, downsampled 0-30-60-90 min / 4 knots) run only if their
#      model files are missing or REFIT=1. The reference fit takes ~1-2 h; to run the fits as
#      separate jobs: RUN=K562_trunc sbatch scripts/slurm/fit.sbatch (and RUN=K562_ds).
#   02 curves, observed means, onset/peak (ChronoPack)   03 figure (R)   04 comparison with v1
#
# Run on a compute node: sbatch scripts/slurm/fig3_v2_pipeline.sbatch (from this directory).
set -euo pipefail

cd "$(dirname "$0")"

CP_ENV="${CP_ENV:-/new-stg/home/alex/miniconda3/envs/chronopack}"
R_ENV="${R_ENV:-/new-stg/home/alex/miniconda3/envs/tradeseq_ce}"
CP_PYTHON="$CP_ENV/bin/python"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/output/.mplconfig}"

step() { echo; echo "=== $* ($(date '+%F %T')) ==="; }

step "00 broad/fine dynamic-mode examples (r9; Fig3_DynamicModes step 01, unchanged)"
FIGURE_VERSION=v4 "$R_ENV/bin/Rscript" scripts/00_plot_dynamic_modes_broad_fine_examples.R K562a

for spec in K562_trunc:8 K562_ds:4; do
  run=${spec%%:*}; k=${spec##*:}
  if [[ "${REFIT:-0}" == "1" || ! -f output/chronopack/$run/sce.1trajectory.nknots$k.rds ]]; then
    step "01 ChronoPack fit $run"
    PATH="$CP_ENV/bin:$PATH" "$CP_PYTHON" scripts/01_fit_downsampling.py --run "$run"
  else
    echo "Using existing fit output/chronopack/$run/sce.1trajectory.nknots$k.rds (REFIT=1 to refit)"
  fi
done

step "02 curves, observed means, onset/peak"
PATH="$CP_ENV/bin:$PATH" "$CP_PYTHON" scripts/02_predict_and_characterize.py

step "03 timing composite"
REQUIRE_QVALUE=true USE_BREAK_AXIS=false OUTPUT_SUFFIX_EXTRA=_q005_no_ybreak \
  "$R_ENV/bin/Rscript" scripts/03_make_timing_composite.R

step "04 comparison with v1"
PATH="$CP_ENV/bin:$PATH" "$CP_PYTHON" scripts/04_compare_to_v1.py

echo "FIG3_V2_PIPELINE_OK"
