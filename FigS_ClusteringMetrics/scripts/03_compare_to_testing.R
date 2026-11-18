#!/usr/bin/env Rscript
# Compare the r9 metrics (step 01) with TESTING.ipynb's saved results (input/testing/, written 2026-09-28 by the
# notebook from the Sep 19 verification run, whose scale factors included the K562_puro rows).
# Reports, per direction and K, the TESTING value, the r9 value and the difference for mean silhouette,
# median / q10 ARI, median mean-cluster Jaccard and per-cluster median Jaccard. Writes output/comparison_to_testing.txt.
args <- commandArgs(trailingOnly = FALSE)
source(file.path(dirname(sub("^--file=", "", args[grep("^--file=", args)])), "common.R"))

new_ce <- readRDS(file.path(RDS_DIR, "clust_eval_ls.rds"))
new_st <- readRDS(file.path(RDS_DIR, "stab_ls.rds"))
old_ce <- readRDS(file.path(TESTING_DIR, "cluster_eval_ls.2rds"))
old_st <- readRDS(file.path(TESTING_DIR, "stab_ls.2rds"))

lines <- c("TESTING.ipynb saved results (input/testing/*.2rds, 2026-09-28) vs r9 recompute (output/rds/)", "")
add <- function(title, df) {
  lines <<- c(lines, title, capture.output(print(as.data.frame(df), row.names = FALSE, digits = 4)), "")
}
dirs <- intersect(names(old_ce), names(new_ce))

sil <- bind_rows(lapply(dirs, function(n) {
  inner_join(old_ce[[n]]$metrics %>% select(K, TESTING = mean_silhouette),
             new_ce[[n]]$metrics %>% select(K, r9 = mean_silhouette), by = "K") %>%
    mutate(direction = n, .before = 1)
})) %>% filter(K <= 20) %>% mutate(diff = r9 - TESTING)
add("Mean silhouette (K <= 20)", sil)

old_k <- sort(unique(old_st[[dirs[1]]]$summary$K))
lines <- c(lines, paste0("TESTING stability K grid: ", paste(old_k, collapse = ", "),
                         "; r9 grid: ", min(K_GRID_STAB), "-", max(K_GRID_STAB),
                         " (same seed and sampling, so K values in both grids are directly comparable)"), "")
stab <- bind_rows(lapply(dirs, function(n) {
  inner_join(old_st[[n]]$summary %>% select(K, TESTING_median_ARI = median_ARI, TESTING_q10_ARI = q10_ARI,
                                            TESTING_median_Jaccard = median_Jaccard),
             new_st[[n]]$summary %>% select(K, r9_median_ARI = median_ARI, r9_q10_ARI = q10_ARI,
                                            r9_median_Jaccard = median_Jaccard), by = "K") %>%
    mutate(direction = n, .before = 1)
}))
add("Bootstrap ARI and mean cluster Jaccard", stab)

if (all(vapply(dirs, function(n) !is.null(old_st[[n]]$cluster_summary), logical(1)))) {
  cj <- bind_rows(lapply(dirs, function(n) {
    inner_join(old_st[[n]]$cluster_summary %>% select(K, cluster, TESTING = median_Jaccard),
               new_st[[n]]$cluster_summary %>% select(K, cluster, r9 = median_Jaccard), by = c("K", "cluster")) %>%
      mutate(direction = n, .before = 1)
  })) %>% filter(K %in% K_CLUSTER_PANEL) %>% mutate(diff = r9 - TESTING, cluster = as.integer(cluster)) %>%
    arrange(direction, K, cluster)
  add("Per-cluster median Jaccard (K = 2-6; cluster ids are each run's own cutree labels)", cj)
} else {
  lines <- c(lines, "TESTING stab_ls has no cluster_summary; per-cluster Jaccard not compared.", "")
}

writeLines(lines, file.path(OUTPUT_DIR, "comparison_to_testing.txt"))
cat(lines, sep = "\n")
cat("COMPARE_TESTING_OK\n")
