#!/usr/bin/env Rscript
# TESTING.ipynb's two clustering-metric figures (execution counts 538 and 539), drawn from step 01's r9 results:
#   clustering_metrics_by_K:  mean silhouette (K <= 20) | bootstrap ARI | mean cluster Jaccard, induced vs repressed
#   cluster_jaccard_by_K:     per-cluster Jaccard (median, 10th percentile) for K = 2-6
# Plot code as in the notebook (except a single shared legend in the first figure); cluster_jac_df (undefined in the saved notebook) is the stability cluster_summary.
args <- commandArgs(trailingOnly = FALSE)
source(file.path(dirname(sub("^--file=", "", args[grep("^--file=", args)])), "common.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

clust_eval.ls <- readRDS(file.path(RDS_DIR, "clust_eval_ls.rds"))
stab.ls <- readRDS(file.path(RDS_DIR, "stab_ls.rds"))
stab_induced <- stab.ls[["induced"]]
stab_repressed <- stab.ls[["repressed"]]

save_both <- function(p, stem, width, height) {
  for (ext in c("pdf", "png")) {
    f <- file.path(FIGURE_DIR, paste0(stem, ".", ext))
    ggsave(f, p, width = width, height = height, dpi = 200)
    message("wrote ", f)
  }
}

# ---- execution count 538 ----
sil_df <- bind_rows(
  induced = clust_eval.ls[["induced"]]$metrics,
  repressed = clust_eval.ls[["repressed"]]$metrics,
  .id = "direction"
) %>% filter(K <= 20)

ari_df <- bind_rows(
  induced = stab_induced$summary,
  repressed = stab_repressed$summary,
  .id = "direction"
)

jac_df <- bind_rows(
  induced = stab_induced$raw,
  repressed = stab_repressed$raw,
  .id = "direction"
) %>%
  group_by(direction, K) %>%
  summarise(
    median_Jaccard = median(mean_Jaccard, na.rm = TRUE),
    q10_Jaccard = quantile(mean_Jaccard, 0.1, na.rm = TRUE),
    .groups = "drop"
  )

base_theme <- theme_bw() +
  theme(legend.position = "right")

p_sil <- ggplot(sil_df, aes(K, mean_silhouette, colour = direction)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 2:20) +
  labs(x = "K", y = "Mean silhouette", colour = NULL) +
  base_theme

p_ari <- ggplot(ari_df, aes(K, median_ARI, colour = direction)) +
  geom_line() + geom_point() +
  geom_line(aes(y = q10_ARI), linetype = 2) +
  scale_x_continuous(breaks = 2:20) +
  labs(x = "K", y = "ARI", colour = NULL,
       subtitle = "Solid = median; dashed = 10th percentile") +
  base_theme

p_jac <- ggplot(jac_df, aes(K, median_Jaccard, colour = direction)) +
  geom_line() + geom_point() +
  geom_line(aes(y = q10_Jaccard), linetype = 2) +
  scale_x_continuous(breaks = 2:20) +
  labs(x = "K", y = "Mean cluster Jaccard", colour = NULL,
       subtitle = "Solid = median; dashed = 10th percentile") +
  base_theme

# one shared legend (the notebook's guides = "collect" kept two copies, from the silhouette and the stability panels)
p538 <- (p_sil | (p_ari + guides(colour = "none")) | (p_jac + guides(colour = "none"))) +
  patchwork::plot_layout(guides = "collect")
save_both(p538, "clustering_metrics_by_K", 16, 5)

# ---- execution count 539 ----
cluster_jac_df <- bind_rows(
  induced = stab_induced$cluster_summary,
  repressed = stab_repressed$cluster_summary,
  .id = "direction"
) %>%
  mutate(cluster = factor(cluster, levels = as.character(sort(unique(as.integer(cluster))))))

p539 <- ggplot(cluster_jac_df %>% filter(K %in% K_CLUSTER_PANEL),
       aes(cluster, median_Jaccard, colour = direction)) +
  geom_point(size = 2.5) +
  geom_segment(aes(xend = cluster, y = q10_Jaccard, yend = median_Jaccard),
               linewidth = 0.7) +
  facet_wrap(~ K, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Reference cluster",
    y = "Jaccard",
    colour = NULL,
    subtitle = "Point = median; lower endpoint = 10th percentile"
  ) +
  theme_bw()
save_both(p539, "cluster_jaccard_by_K", 16, 5)

cat("PLOT_METRICS_OK\n")
