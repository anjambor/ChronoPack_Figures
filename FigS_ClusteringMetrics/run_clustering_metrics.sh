#!/usr/bin/env bash
# FigS_ClusteringMetrics: gene-trajectory clustering metrics (silhouette, bootstrap ARI / Jaccard) on ChronoPack r9.
# Run on a compute node:  sbatch scripts/slurm/clustering_metrics.sbatch
set -euo pipefail
cd "$(dirname "$0")"
RSCRIPT="${RSCRIPT_BIN:-$HOME/miniconda3/envs/chronopack/bin/Rscript}"
mkdir -p output/logs
"$RSCRIPT" scripts/01_clustering_metrics.R
"$RSCRIPT" scripts/02_plot_metrics.R
"$RSCRIPT" scripts/03_compare_to_testing.R
"$RSCRIPT" scripts/04_plot_composite.R
echo CLUSTERING_METRICS_PIPELINE_OK
