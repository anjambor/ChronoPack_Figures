# FigS_ClusteringMetrics: gene-trajectory clustering metrics on ChronoPack r9

This folder recomputes the clustering-metric panels of
`chrono_seq/analyses/k562-tnfa-puromycin.hg38/notebooks/TESTING.ipynb` (execution counts 538 and 539) from the
**ChronoPack r9 gene clustering**:
- mean silhouette;
- bootstrap ARI;
- mean cluster Jaccard;
- per-cluster Jaccard.

Each metric is computed for the induced and repressed gene-trajectory clusterings. Nothing in `chrono_seq/` or
`ChronoPack/` is changed.

**Output:**
- `output/figures/clustering_metrics_by_K.{pdf,png}`: mean silhouette (K ≤ 20) | ARI | mean cluster Jaccard. This is
  execution count 538.
- `output/figures/cluster_jaccard_by_K.{pdf,png}`: per-cluster Jaccard, K = 2–6. This is execution count 539.
- **`output/figures/clustering_metrics_composite.{pdf,png}`**: the publication composite, 180 × 130 mm, combining both
  figures (step 04; layout and caption below).

### Composite layout

- **Top row, (a)–(c):** mean silhouette, bootstrap ARI and mean cluster Jaccard against K. The three panels share the
  K axis.
  - (b) and (c) share a 0.3–1 scale, since both are agreement scores.
  - Dotted lines mark the clusterings r9 uses: broad K = 3 (grey) and fine K = the elbow K (induced 17, repressed 16).
- **Bottom row, (d):** per-cluster Jaccard for K = 2–6, laid out as direction rows × K columns, with column widths
  proportional to K.
  - Induced and repressed are separate trees, so their cluster numbers do not correspond. That is why they get
    separate rows; the notebook's version overlaid them on one axis.
- **Legend:** one shared legend (colour = direction, linetype = statistic) replaces the notebook's subtitles.

**Figure legend (draft):**

> **Figure S#. Stability of gene-trajectory clustering across the number of clusters.**
> Significant induced (red; n = 2,424) and repressed (blue; n = 2,667) gene trajectories, pooled across the K562a,
> K562b and THP1 datasets, were clustered separately by Ward's hierarchical clustering. The clustering used row-z-scored
> fitted trajectories together with their time derivatives, with the derivatives weighted α = 0.1.
> **(a)** Mean silhouette width of the full-data clustering cut at K = 2–20 clusters.
> **(b, c)** Agreement between the full-data clustering and clusterings of 200 random subsamples, each containing 80% of
> the trajectories. Each subsample was compared with the full-data partition restricted to the same trajectories, by the
> adjusted Rand index (b) and by the mean over clusters of each full-data cluster's best-matching Jaccard index (c).
> Solid lines show the median over subsamples and dashed lines show the 10th percentile.
> In (a)–(c), dotted vertical lines mark the broad partition (K = 3; grey) and the fine partitions chosen by the elbow
> criterion (induced K = 17, repressed K = 16; coloured).
> **(d)** Jaccard index of each full-data cluster for K = 2–6. Points show the median over the 200 subsamples, and lines
> extend down to the 10th percentile. The dashed line marks Jaccard = 0.5. Cluster numbers are labels from the
> full-data tree and do not correspond between induced and repressed trajectories.


## Why rerun

TESTING.ipynb's panels do not come from r9's inputs:
- **Deleted inputs.** The notebook reads `X`, `Xd` and `pc_df` from
  `~/tmp/ChronoPack/verification_runs/manuscript_full_pipeline_20260919_165539/`. That Sep 19 run has since been deleted.
- **Scale factors included K562_puro.** That run had 6,723 trajectories, including K562_puro. The notebook drops
  K562_puro only when it picks the induced and repressed genes; the scale factors `sx` and `sdx` were computed over all
  6,723 rows. r9 has no K562_puro, so its scale factors use only its 5,091 trajectories.
  - This shifts the relative weight of the derivative block by about 0.5% (estimated with r8, which also has K562_puro).
  - The notebook's saved results show the effect is small (see Results). Mean silhouette changes by about 0.0003 at
    most K, and by 0.019 at repressed K = 11.
- **Silhouette came from an earlier run.** The silhouette panel used a `clust_eval.ls` built by an earlier run of the
  clustering cells. Those cells were rerun later, so its settings are not saved.
- **Per-cluster panel input is undefined.** `cluster_jac_df`, the per-cluster panel's input, is not defined anywhere in
  the saved notebook.

## Method

**Clustering matrix.** It is built exactly as `cp.tl.cluster_gene_trajectories` builds it (the package R helper
`chronopack_cluster_trajectory_matrix_to_files`), from r9's `gene_trajectory_fit.csv` and time grid:
- **Reference:** the significant trajectories, pvalue_1 < 0.05. That is all 5,091 r9 rows: K562a 2,022, K562b 1,950 and
  THP1 1,119.
- **Matrices:** X is the per-gene z-scored curves. Xd is the per-gene z-scored finite-difference derivative.
- **Scaling:** each matrix is divided by `sqrt(sum of column variances)` over the reference.
- **Joint matrix:** `cbind(sqrt(1 - α) X_n, sqrt(α) Xd_n)` with α = 0.1.
- **Subsets:** induced = peak_log2fc > 0 (2,424 genes); repressed = peak_log2fc < 0 (2,667 genes).

**Metrics:**
- **Mean silhouette:** Ward.D2 clustering on Euclidean distances, then the package's
  `calculate_trajectory_clustering_metrics` over K = 2–100, as in r9. The figure shows K ≤ 20.
- **ARI, mean cluster Jaccard, per-cluster Jaccard:** TESTING.ipynb's `cluster_stability`, copied verbatim into
  `scripts/common.R`.
  - 200 bootstrap repeats; each re-clusters a random 80% of genes. The seed is 1.
  - Each repeat's partition is compared with the full-data tree at the same K.
  - Mean cluster Jaccard is the mean, over reference clusters, of each cluster's best Jaccard match in the repeat.
  - Results are summarised as the median and 10th percentile over the 200 repeats.
  - The notebook ran K = 2–6; this folder runs K = 2–20. The sampling does not depend on the K grid, so K = 2–6 are
    exactly what a K = 2–6 run would give.

**Checks** (step 01; written to `output/tables/input_checks.txt`; the run fails if any check fails):
- The z-scored matrices equal r9's `gene_trajectory_fit_{,ddt_}significant_z.csv`.
- The silhouettes equal r9's `cluster_metrics.{induced,repressed}.csv`.
- The K = 3 (broad) and elbow-K (fine) partitions equal r9's `cluster_{broad,fine}_{induced,repressed}` in
  `pc_df.new.csv`.

## Results (2026-10-21)

All step 01 checks passed:
- The z-scored matrices match r9's to within 6e-11.
- The silhouettes match r9's to within 5e-16.
- The broad and fine partitions match r9's: induced k = 3 and 17, repressed k = 3 and 16.

The stability bootstrap took 7 min for the induced genes and 10 min for the repressed genes.

**r9 values** (`output/tables/`):

| | K = 2 | K = 3 | K = 4 | K = 5 | K = 6 |
|---|---|---|---|---|---|
| Induced mean silhouette | 0.564 | 0.455 | 0.456 | 0.460 | 0.359 |
| Induced median ARI | 0.84 | 0.89 | 0.88 | 0.77 | 0.59 |
| Induced median mean Jaccard | 0.90 | 0.87 | 0.74 | 0.71 | 0.62 |
| Repressed mean silhouette | 0.447 | 0.512 | 0.465 | 0.458 | 0.438 |
| Repressed median ARI | 0.80 | 0.75 | 0.62 | 0.60 | 0.66 |
| Repressed median mean Jaccard | 0.89 | 0.80 | 0.58 | 0.61 | 0.61 |

**Compared with TESTING.ipynb's saved results** (`output/comparison_to_testing.txt`; `stab_ls.2rds` /
`cluster_eval_ls.2rds`, 2026-09-28):
- **Mean silhouette (K ≤ 20):** within 0.0003 everywhere except repressed K = 11, where the notebook has 0.427 and r9
  has 0.447. At that K, a near-tied Ward merge falls the other way.
- **ARI and mean cluster Jaccard (K = 2–6, the notebook's grid):** within 0.04. The largest gaps are repressed K = 2
  (ARI 0.84 vs 0.80) and induced K = 2 (0.86 vs 0.84). Small input differences change which trees the bootstrap
  samples produce, so these values move slightly even though the seed is the same.
- **Per-cluster median Jaccard (K = 2–6):** within 0.025, except three repressed K = 6 clusters (up to 0.054 lower in
  r9).

The shape of every curve is unchanged:
- **Induced:** the mean silhouette peaks at K = 2. The median ARI stays at 0.84–0.89 for K = 2–4, then falls to 0.77 at
  K = 5 and 0.59 at K = 6.
- **Repressed:** the mean silhouette peaks at K = 3. The median ARI is 0.80 and 0.75 at K = 2 and 3, then drops to 0.62
  at K = 4.

## Run

From this directory:

```bash
sbatch scripts/slurm/clustering_metrics.sbatch
```

Logs go to `output/logs/`. A successful run ends with `CLUSTERING_METRICS_PIPELINE_OK`.

| Step | Script | What it does |
|---|---|---|
| 01 | `01_clustering_metrics.R` | Builds the matrix and runs the checks, the hierarchical metrics and the bootstrap stability. Writes `output/tables/{silhouette,stability_summary,stability_raw,cluster_jaccard_summary,cluster_jaccard_raw}.csv` and `output/rds/{clust_eval_ls,stab_ls}.rds` (the notebook's object structure). |
| 02 | `02_plot_metrics.R` | The notebook's two figures. The plot code is the notebook's, except that the first figure has one shared legend instead of two. |
| 03 | `03_compare_to_testing.R` | Compares the results with TESTING.ipynb's saved results and writes `output/comparison_to_testing.txt`. |
| 04 | `04_plot_composite.R` | The publication composite. It checks that its medians equal step 01's `stability_summary.csv`. |

`scripts/common.R` holds the paths, the constants and the notebook's stability functions.

## Inputs (`input/`)

- `r9_gene_clustering/` is a symlink to `ChronoPack/verification_runs/full_pipeline_analyses_20260930_r9/data/gene_clustering`.
  It provides `gene_trajectory_fit.csv`, `gene_trajectory_time_grid.csv` and `pc_df.new.csv`, plus the files the checks
  read.
- `testing/` holds `stab_ls.2rds` and `cluster_eval_ls.2rds`, copied from `chrono_seq/.../notebooks/` (written by
  TESTING.ipynb on 2026-09-28). Only step 03 reads them.

## External dependencies

- `~/miniconda3/envs/chronopack` Rscript (dplyr, tidyr, ggplot2, patchwork, matrixStats). You can override it with
  `RSCRIPT_BIN`.
- `../ChronoPack/src/chronopack/_backends/r/dimred_clustering_helpers.R`, the package's R helpers. They are sourced
  read-only.
