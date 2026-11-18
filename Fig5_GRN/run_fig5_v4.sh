#!/usr/bin/env bash
# Fig5_GRN_v4 pipeline: Fig5_GRN_v3 (BEELINE comparators from the r9 single-cell P < 0.05 TFs + 2000 run, motif-prior
# methods in a second method set) with ChronoPack's TF rankings computed by the current ChronoPack (rank_regulators
# defaults, commit e0c89fb) from the pinned OmniPath snapshot.
# ChronoPack GRN inference (r9, once) + benchmark tables + publication-v4 figures for each method set:
#   expression  (FIG5_METHOD_SET default) -> output/benchmark_outputs/
#   with_prior                             -> output/with_prior/benchmark_outputs/
#
#   LIVE_OMNIPATH=1   rank TFs from the current OmniPath (network access) instead of the snapshot.
#   FIG5_OUTPUT_DIR   redirect all outputs (default: ./output).
#
# Run on a compute node: sbatch scripts/slurm/fig5_pipeline.sbatch (from this directory).
set -euo pipefail

cd "$(dirname "$0")"

CP_ENV="${CP_ENV:-/new-stg/home/alex/miniconda3/envs/chronopack}"
CP_PYTHON="${CP_PYTHON:-$CP_ENV/bin/python}"
export PATH="$CP_ENV/bin:$PATH"  # R for rpy2 (score_targets)
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/output/.mplconfig}"

step() { echo; echo "=== $* ($(date '+%F %T')) ==="; }

step "01 ChronoPack GRN (cp.tl.grn)"
if [[ "${LIVE_OMNIPATH:-0}" == "1" ]]; then
  "$CP_PYTHON" scripts/01_run_chronopack_grn.py --live-omnipath
else
  "$CP_PYTHON" scripts/01_run_chronopack_grn.py
fi

for set in expression with_prior; do
  export FIG5_METHOD_SET=$set
  step "02 benchmark tables + figures [$set]"
  "$CP_PYTHON" scripts/02_reproduce_figures.py

  step "02b 5d expression curves from the ChronoPack 1-lineage K562 fit [$set]"
  "$CP_PYTHON" scripts/02b_plot_expr_curves_1lineage.py

  step "02c 5b K562 ChIP strip-forests with EPR at k = 500 [$set]"
  "$CP_PYTHON" scripts/02c_k562_chip_epr_at_500.py
done

# with_prior on the shared gene universe without Replogle / Yao -> with_prior figures 5{b,c,d,e}_v2 (after the
# full with_prior run, which clears output/with_prior/benchmark_outputs/)
export FIG5_METHOD_SET=with_prior FIG5_UNIVERSE=shared
for s in 02_reproduce_figures 02b_plot_expr_curves_1lineage 02c_k562_chip_epr_at_500 02d_install_shared_universe_v2; do
  step "$s [with_prior, shared universe]"
  "$CP_PYTHON" "scripts/$s.py"
done
unset FIG5_UNIVERSE

step "02e Perturb-seq (Replogle, Yao) as functional references [with_prior]"
"$CP_PYTHON" scripts/02e_perturbseq_reference.py

export FIG5_METHOD_SET=expression
# Step 03 (checks against Fig5_GRN_v3) was removed 2026-11-16 (v3 is in ../old/; user request: no dependencies on
# old/). Its last report, output/comparison_to_v3.txt (2026-10-22, all checks passed), is kept.

step "05 composite figure (layout of ref/FigureComposite.v6.png, expression set)"
"$CP_PYTHON" scripts/05_composite_figure.py

echo "FIG5_V4_PIPELINE_OK"
