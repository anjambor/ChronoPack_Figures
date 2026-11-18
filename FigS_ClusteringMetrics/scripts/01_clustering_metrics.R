#!/usr/bin/env Rscript
# Clustering metrics for the induced / repressed gene-trajectory clustering on ChronoPack r9 inputs.
#
# Clustering matrix: built as cp.tl.cluster_gene_trajectories builds it (chronopack_cluster_trajectory_matrix_to_files
# in the package's dimred_clustering_helpers.R), from r9's gene_trajectory_fit.csv and time grid:
#   reference = the significant trajectories (pvalue_1 < 0.05; all 5,091 r9 rows, no K562_puro);
#   X = z_row(curves), Xd = z_row(finite-difference derivative); each divided by its total column variance over
#   the reference; X_joint = cbind(sqrt(1 - alpha) X_n, sqrt(alpha) Xd_n), alpha = 0.1;
#   induced = reference rows with peak_log2fc > 0, repressed = peak_log2fc < 0.
#
# 1. Hierarchical metrics: Ward.D2 on Euclidean distances, calculate_trajectory_clustering_metrics (package) over
#    K = 2-100 -> mean silhouette. Checked against r9's cluster_metrics.{induced,repressed}.csv and r9's broad (k = 3)
#    and fine (elbow k) partitions in pc_df.new.csv.
# 2. Bootstrap stability: TESTING.ipynb's cluster_stability (B = 200, 80% of genes, seed 1) over K = 2-20 ->
#    ARI, mean cluster Jaccard and per-cluster Jaccard against the full-data tree.
#
# Writes output/tables/{silhouette,stability_summary,stability_raw,cluster_jaccard_summary,cluster_jaccard_raw}.csv,
# output/tables/input_checks.txt and output/rds/{clust_eval_ls,stab_ls}.rds.
args <- commandArgs(trailingOnly = FALSE)
source(file.path(dirname(sub("^--file=", "", args[grep("^--file=", args)])), "common.R"))
source(CHRONOPACK_R_HELPERS)

checks <- character()
note <- function(...) {
  msg <- paste0(...)
  message(msg)
  checks <<- c(checks, msg)
}

traj_mat <- as.matrix(read.csv(file.path(R9_DIR, "gene_trajectory_fit.csv"), row.names = 1, check.names = FALSE))
tt <- read.csv(file.path(R9_DIR, "gene_trajectory_time_grid.csv"), check.names = FALSE)$time
pc_df <- read.csv(file.path(R9_DIR, "pc_df.new.csv"), row.names = 1, check.names = FALSE)
stopifnot(setequal(rownames(pc_df), rownames(traj_mat)), !anyDuplicated(rownames(pc_df)), length(tt) == ncol(traj_mat))
pc_df <- pc_df[rownames(traj_mat), , drop = FALSE]  # pc_df.new.csv is in cluster order; align to the trajectory rows
note("r9 trajectories: ", nrow(traj_mat), " x ", ncol(traj_mat), "; cell types: ",
     paste(names(table(pc_df$cell_type)), table(pc_df$cell_type), sep = " ", collapse = ", "))
stopifnot(!any(pc_df$cell_type == "K562_puro"))

# ---- clustering matrix, as chronopack_cluster_trajectory_matrix_to_files ----
sig <- rownames(pc_df)[!is.na(pc_df$pvalue_1) & pc_df$pvalue_1 < PVALUE_THRESHOLD]
ref_mat <- traj_mat[sig, , drop = FALSE]
X_ref <- z_row(ref_mat)
Xd_ref <- ref_mat[, -1, drop = FALSE] - ref_mat[, -ncol(ref_mat), drop = FALSE]
Xd_ref <- sweep(Xd_ref, 2, diff(tt), "/")
Xd_ref <- z_row(Xd_ref)
inputs <- prepare_dimred_clustering_inputs(X_ref, Xd_ref)
X_all <- cbind(sqrt(1 - ALPHA) * inputs$X_n, sqrt(ALPHA) * inputs$Xd_n)
note("reference (significant) rows: ", length(sig),
     "; scale factors sx = ", signif(sqrt(sum(apply(X_ref, 2, var))), 6),
     ", sdx = ", signif(sqrt(sum(apply(Xd_ref, 2, var))), 6))

# check: the z-scored matrices equal r9's written ones
for (f in c("gene_trajectory_fit_significant_z.csv", "gene_trajectory_fit_ddt_significant_z.csv")) {
  M <- as.matrix(read.csv(file.path(R9_DIR, f), row.names = 1, check.names = FALSE))
  mine <- if (grepl("ddt", f)) Xd_ref else X_ref
  d <- max(abs(M[rownames(mine), ] - mine))
  note(f, ": max |diff| = ", signif(d, 3))
  stopifnot(identical(sort(rownames(M)), sort(rownames(mine))), d < 1e-8)
}

ids <- list(
  induced = sig[pc_df[sig, "peak_log2fc"] > LOG2FC_THRESHOLD],
  repressed = sig[pc_df[sig, "peak_log2fc"] < -LOG2FC_THRESHOLD]
)

same_partition <- function(a, b) {
  tab <- table(a, b)
  all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1)
}

clust_eval.ls <- list()
stab.ls <- list()
for (name in names(ids)) {
  X <- X_all[ids[[name]], , drop = FALSE]
  note("\n", name, ": ", nrow(X), " genes")

  # 1. hierarchical metrics (as the package)
  D <- dist(X, method = "euclidean")
  hc <- hclust(D, method = "ward.D2")
  ce <- calculate_trajectory_clustering_metrics(X_mat = X, D = D, hc = hc, K_grid = K_GRID_HIER,
                                                wss_fraction = WSS_FRACTION, select_k_by = "elbow")
  clust_eval.ls[[name]] <- ce
  r9m <- read.csv(file.path(R9_DIR, paste0("cluster_metrics.", name, ".csv")))
  d <- max(abs(r9m$mean_silhouette[match(ce$metrics$K, r9m$K)] - ce$metrics$mean_silhouette))
  note("  mean silhouette vs r9 cluster_metrics.", name, ".csv: max |diff| = ", signif(d, 3),
       "; elbow k = ", ce$selection$k_elbow)
  stopifnot(d < 1e-8)
  for (lvl in c("broad", "fine")) {
    k <- if (lvl == "broad") 3L else ce$selection$k_elbow
    r9_lab <- pc_df[ids[[name]], paste0("cluster_", lvl, "_", name)]
    ok <- same_partition(cutree(hc, k = k), r9_lab)
    note("  k = ", k, " partition equals r9 cluster_", lvl, "_", name, ": ", ok)
    stopifnot(ok)
  }

  # 2. bootstrap stability (TESTING.ipynb)
  t0 <- Sys.time()
  stab.ls[[name]] <- cluster_stability(X, K_grid = K_GRID_STAB, B = B, frac = FRAC, seed = SEED)
  note("  stability: B = ", B, ", frac = ", FRAC, ", K = ", min(K_GRID_STAB), "-", max(K_GRID_STAB),
       " (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
}

bind <- function(f) bind_rows(lapply(names(ids), function(n) mutate(f(n), direction = n, .before = 1)))
write.csv(bind(function(n) clust_eval.ls[[n]]$metrics), file.path(TABLE_DIR, "silhouette.csv"), row.names = FALSE)
write.csv(bind(function(n) stab.ls[[n]]$summary), file.path(TABLE_DIR, "stability_summary.csv"), row.names = FALSE)
write.csv(bind(function(n) stab.ls[[n]]$raw), file.path(TABLE_DIR, "stability_raw.csv"), row.names = FALSE)
write.csv(bind(function(n) stab.ls[[n]]$cluster_summary), file.path(TABLE_DIR, "cluster_jaccard_summary.csv"),
          row.names = FALSE)
write.csv(bind(function(n) stab.ls[[n]]$cluster_raw), file.path(TABLE_DIR, "cluster_jaccard_raw.csv"),
          row.names = FALSE)
saveRDS(clust_eval.ls, file.path(RDS_DIR, "clust_eval_ls.rds"))
saveRDS(stab.ls, file.path(RDS_DIR, "stab_ls.rds"))
writeLines(checks, file.path(TABLE_DIR, "input_checks.txt"))

print(as.data.frame(bind(function(n) stab.ls[[n]]$summary)) %>% filter(K <= 8), digits = 3)
cat("CLUSTERING_METRICS_OK\n")
