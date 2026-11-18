# Shared paths, constants and the TESTING.ipynb stability functions for FigS_ClusteringMetrics.
#
# The stability functions (ari, best_jaccard, cluster_stability) are copied verbatim from
# chrono_seq/analyses/k562-tnfa-puromycin.hg38/notebooks/TESTING.ipynb (cell at execution count 530).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

get_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  normalizePath(file.path(dirname(f), ".."))
}
ROOT <- get_root()
INPUT_DIR <- file.path(ROOT, "input")
R9_DIR <- file.path(INPUT_DIR, "r9_gene_clustering")  # -> ChronoPack r9 data/gene_clustering
TESTING_DIR <- file.path(INPUT_DIR, "testing")         # TESTING.ipynb's saved stab_ls / cluster_eval_ls (copied)
OUTPUT_DIR <- file.path(ROOT, "output")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
RDS_DIR <- file.path(OUTPUT_DIR, "rds")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")

# chronopack package R helpers (z_row, scale_by_total_column_variance, calculate_trajectory_clustering_metrics)
CHRONOPACK_R_HELPERS <- normalizePath(file.path(
  ROOT, "..", "..", "ChronoPack", "src", "chronopack", "_backends", "r", "dimred_clustering_helpers.R"))

# As in r9's 01_DimensionalityReduction&Clustering notebook (cp.tl.cluster_gene_trajectories)
ALPHA <- 0.1
PVALUE_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 0
K_GRID_HIER <- 2:100   # r9's grid; the silhouette panel shows K <= 20, as in TESTING.ipynb
WSS_FRACTION <- 0.95

# As in TESTING.ipynb's stability cells (B, frac, seed), over K = 2-20
K_GRID_STAB <- 2:20
B <- 200
FRAC <- 0.8
SEED <- 1
K_CLUSTER_PANEL <- 2:6  # per-cluster Jaccard panel (TESTING.ipynb execution count 539)

for (d in c(TABLE_DIR, RDS_DIR, FIGURE_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- TESTING.ipynb stability functions (verbatim) ----

ari <- function(a, b) {
  z <- table(a, b)
  c2 <- function(x) x * (x - 1) / 2
  nij <- sum(c2(z))
  ai <- sum(c2(rowSums(z)))
  bj <- sum(c2(colSums(z)))
  ex <- ai * bj / c2(sum(z))
  den <- 0.5 * (ai + bj) - ex
  if (den == 0) 1 else (nij - ex) / den
}

best_jaccard <- function(ref, test) {
  cls <- sort(unique(ref))
  setNames(sapply(cls, function(cl) {
    A <- ref == cl
    max(sapply(unique(test), function(cl2) {
      B <- test == cl2
      sum(A & B) / sum(A | B)
    }))
  }), cls)
}

cluster_stability <- function(X, K_grid = 2:20, B = 200,
                              frac = 0.8, seed = 1) {
  X <- as.matrix(X)
  set.seed(seed)

  hc0 <- hclust(dist(X), method = "ward.D2")
  ref <- cutree(hc0, k = K_grid)
  if (is.null(dim(ref))) ref <- matrix(ref, ncol = 1)

  out <- vector("list", B)
  cluster_out <- vector("list", B)

  for (b in seq_len(B)) {
    ii <- sort(sample(seq_len(nrow(X)), ceiling(frac * nrow(X))))
    hcb <- hclust(dist(X[ii, , drop = FALSE]), method = "ward.D2")
    clb <- cutree(hcb, k = K_grid)
    if (is.null(dim(clb))) clb <- matrix(clb, ncol = 1)

    tmp <- lapply(seq_along(K_grid), function(j) {
      r <- ref[ii, j]
      q <- clb[, j]
      jac <- best_jaccard(r, q)

      list(
        global = data.frame(
          replicate = b,
          K = K_grid[j],
          ARI = ari(r, q),
          mean_Jaccard = mean(jac),
          min_Jaccard = min(jac)
        ),
        cluster = data.frame(
          replicate = b,
          K = K_grid[j],
          cluster = names(jac),
          Jaccard = as.numeric(jac)
        )
      )
    })

    out[[b]] <- do.call(rbind, lapply(tmp, `[[`, "global"))
    cluster_out[[b]] <- do.call(rbind, lapply(tmp, `[[`, "cluster"))
  }

  raw <- do.call(rbind, out)
  cluster_raw <- do.call(rbind, cluster_out)

  summary <- raw %>%
    group_by(K) %>%
    summarise(
      median_ARI = median(ARI),
      q10_ARI = quantile(ARI, 0.10),
      median_Jaccard = median(mean_Jaccard),
      median_min_Jaccard = median(min_Jaccard),
      .groups = "drop"
    )

  cluster_summary <- cluster_raw %>%
    group_by(K, cluster) %>%
    summarise(
      median_Jaccard = median(Jaccard),
      q10_Jaccard = quantile(Jaccard, 0.10),
      .groups = "drop"
    )

  list(
    raw = raw,
    summary = summary,
    cluster_raw = cluster_raw,
    cluster_summary = cluster_summary
  )
}
