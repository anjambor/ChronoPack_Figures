#!/usr/bin/env bash
# Submit from FigS_ClusteringMetrics/:  sbatch scripts/slurm/clustering_metrics.sbatch
#SBATCH --job-name=figS_clustering_metrics
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the FigS_ClusteringMetrics directory}"
[[ -f run_clustering_metrics.sh ]] || { echo "Submit from the FigS_ClusteringMetrics directory" >&2; exit 2; }
./run_clustering_metrics.sh
