#!/usr/bin/env bash
# One method on one dataset. Submit from Run_PriorGRN/ (run_prior_grn.sh sets the resources per method):
#   METHOD=PYSCENIC|CELLORACLE|INFERELATOR DATASET=K562-r9-sc-cpfilt-p05-TFs2000 sbatch scripts/slurm/method.sbatch
#SBATCH --job-name=priorgrn
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Run_PriorGRN directory}"
[[ -f run_prior_grn.sh ]] || { echo "Submit from the Run_PriorGRN directory" >&2; exit 2; }
: "${METHOD:?}" "${DATASET:?}"
N=${SLURM_CPUS_PER_TASK:-16}
IN=output/inputs/$DATASET
OUT=output/networks/$DATASET/$METHOD
mkdir -p "$OUT"
export PYTHONNOUSERSITE=1
SIF="singularity exec --cleanenv --env PYTHONNOUSERSITE=1 -B $PWD"
echo "$METHOD $DATASET on $(hostname), $N cpus, start $(date)"
start=$(date +%s)
case $METHOD in
  PYSCENIC)
    # GRNBoost2 without dask (pySCENIC's multiprocessing runner), then cisTarget pruning (pySCENIC protocol options)
    $SIF resources/containers/pyscenic_0.12.1.sif arboreto_with_multiprocessing.py \
      "$IN/expr_log1p.tsv" "$IN/tfs_pyscenic.txt" --method grnboost2 --output "$OUT/adjacencies.tsv" \
      --num_workers "$N" --seed 777
    $SIF resources/containers/pyscenic_0.12.1.sif pyscenic ctx "$OUT/adjacencies.tsv" \
      resources/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
      resources/cistarget/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
      --annotations_fname resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
      --expression_mtx_fname "$IN/expr_log1p.tsv" --output "$OUT/regulons.csv" --mask_dropouts --num_workers "$N"
    $SIF resources/containers/pyscenic_0.12.1.sif python scripts/02_pyscenic_edges.py "$DATASET" ;;
  CELLORACLE)
    $SIF resources/containers/celloracle_0.18.0.sif python scripts/03_celloracle.py "$DATASET" "$N" ;;
  INFERELATOR)
    .conda_env/bin/python scripts/04_inferelator.py "$DATASET" "$N" ;;
  *) echo "unknown METHOD $METHOD" >&2; exit 2 ;;
esac
echo $(( $(date +%s) - start )) > "$OUT/time_seconds.txt"
echo "PRIORGRN_METHOD_OK $METHOD $DATASET $(cat "$OUT/time_seconds.txt") s"
