#!/usr/bin/env bash
# TF-balanced global edge recovery (side analysis; does not change the Fig5 pipeline outputs).
# Submit from Fig5_GRN/:  sbatch scripts/slurm/tf_balanced_global.sbatch
#   STEPS="04a 04" (default) or STEPS=04 to reuse existing DE scores; STEPS=04b for per-TF peak-DE scores;
#   04c / 04d / 04e for the calibration, score-combination and gene-universe side analyses;
#   04f two-tier completion (comparators completed with ChronoPack's candidate genes first);
#   04g gene-set target enrichment and a toy model of random completion.
#SBATCH --job-name=fig5v4_tfbal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Fig5_GRN directory}"
[[ -f run_fig5_v4.sh ]] || { echo "Submit from the Fig5_GRN directory" >&2; exit 2; }
CP_ENV="${CP_ENV:-/new-stg/home/alex/miniconda3/envs/chronopack}"
export PATH="$CP_ENV/bin:$PATH"  # R for rpy2
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/output/.mplconfig}"
for step in ${STEPS:-04a 04}; do
  case "$step" in
    04a) "$CP_ENV/bin/python" scripts/04a_de_target_scores.py ;;
    04)  "$CP_ENV/bin/python" scripts/04_tf_balanced_global_recall.py ;;
    04b) "$CP_ENV/bin/python" scripts/04b_per_tf_de_peak.py ;;
    04c) "$CP_ENV/bin/python" scripts/04c_de_calibration_peak_robustness.py ;;
    04d) "$CP_ENV/bin/python" scripts/04d_score_combination_sensitivity.py ;;
    04e) "$CP_ENV/bin/python" scripts/04e_gene_universe_sensitivity.py ;;
    04f) "$CP_ENV/bin/python" scripts/04f_two_tier_completion.py ;;
    04g) "$CP_ENV/bin/python" scripts/04g_gene_set_enrichment.py ;;
    *) echo "unknown step $step" >&2; exit 2 ;;
  esac
done
echo "TF_BALANCED_OK"
