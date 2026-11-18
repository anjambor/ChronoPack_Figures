#!/usr/bin/env bash
# Submit from Fig5_GRN/:  sbatch scripts/slurm/fig5_shared_universe.sbatch
# with_prior figures on the shared gene universe without Replogle / Yao (FIG5_UNIVERSE=shared) -> 5{b,c,d,e}_v2.
# Needs step 01's output/chronopack/ (shared universe) and is also run at the end of run_fig5_v4.sh.
#SBATCH --job-name=fig5v4_shared
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Fig5_GRN directory}"
[[ -f run_fig5_v4.sh ]] || { echo "Submit from the Fig5_GRN directory" >&2; exit 2; }
CP_ENV="${CP_ENV:-/new-stg/home/alex/miniconda3/envs/chronopack}"
export PATH="$CP_ENV/bin:$PATH"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/output/.mplconfig}"
export FIG5_METHOD_SET=with_prior FIG5_UNIVERSE=shared
for s in 02_reproduce_figures 02b_plot_expr_curves_1lineage 02c_k562_chip_epr_at_500 02d_install_shared_universe_v2; do
  echo; echo "=== $s [with_prior, shared universe] ($(date '+%F %T')) ==="
  python "scripts/$s.py"
done

# Perturb-seq (Replogle K562, Yao THP1) as functional references for the with_prior networks -> perturbseq_reference_v2
echo; echo "=== 02e_perturbseq_reference [with_prior] ($(date '+%F %T')) ==="
unset FIG5_UNIVERSE
python scripts/02e_perturbseq_reference.py
