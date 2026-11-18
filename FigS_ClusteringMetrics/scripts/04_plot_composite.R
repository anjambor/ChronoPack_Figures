#!/usr/bin/env Rscript
# Publication composite of the clustering metrics (step 01's r9 results), 180 mm wide:
#   a  mean silhouette vs K            b  bootstrap ARI vs K            c  mean cluster Jaccard vs K
#   d  per-cluster Jaccard, K = 2-6: rows = direction, columns = K (column width proportional to K)
# a-c share the K axis and mark the clusterings r9 uses (broad K = 3 for both directions; fine K = the elbow K of each
# direction, read from step 01's metrics). b and c share a 0-1 agreement scale. Median = solid, 10th percentile =
# dashed (b, c) or the lower end of each segment (d). In d, cluster numbers are cutree labels of the full-data tree;
# induced and repressed clusters are unrelated, hence separate rows.
# Writes output/figures/clustering_metrics_composite.{pdf,png}.
args <- commandArgs(trailingOnly = FALSE)
source(file.path(dirname(sub("^--file=", "", args[grep("^--file=", args)])), "common.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

clust_eval.ls <- readRDS(file.path(RDS_DIR, "clust_eval_ls.rds"))
stab.ls <- readRDS(file.path(RDS_DIR, "stab_ls.rds"))
DIRS <- c("induced", "repressed")
COLS <- c(induced = "#C0392B", repressed = "#2166AC")
LABS <- c(induced = "Induced", repressed = "Repressed")
K_MAX <- max(K_GRID_STAB)
BROAD_K <- 3

sil_df <- bind_rows(lapply(DIRS, function(n) mutate(clust_eval.ls[[n]]$metrics, direction = n))) %>% filter(K <= K_MAX)
fine_k <- tibble(direction = DIRS, K = vapply(DIRS, function(n) as.numeric(clust_eval.ls[[n]]$selection$k_elbow), 1))
stab_long <- bind_rows(lapply(DIRS, function(n) mutate(stab.ls[[n]]$raw, direction = n))) %>%
  group_by(direction, K) %>%
  summarise(ARI_median = median(ARI), ARI_q10 = quantile(ARI, 0.1),
            Jaccard_median = median(mean_Jaccard), Jaccard_q10 = quantile(mean_Jaccard, 0.1), .groups = "drop") %>%
  pivot_longer(-c(direction, K), names_to = c("metric", "stat"), names_sep = "_") %>%
  mutate(stat = factor(stat, levels = c("median", "q10"), labels = c("Median", "10th percentile")))
# the summary must equal the notebook-style tables written by step 01
chk <- read.csv(file.path(TABLE_DIR, "stability_summary.csv"))
m <- stab_long %>% filter(metric == "ARI", stat == "Median") %>% inner_join(chk, by = c("direction", "K"))
stopifnot(nrow(m) == nrow(chk), max(abs(m$value - m$median_ARI)) < 1e-12)

theme_pub <- theme_classic(base_size = 7) +
  theme(
    axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
    axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 7),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
    plot.tag = element_text(size = 9, face = "bold"),
    legend.title = element_blank(), legend.text = element_text(size = 6.5),
    legend.key.width = unit(5, "mm"), legend.key.height = unit(2.5, "mm"),
    strip.background = element_blank(), strip.text = element_text(size = 6.5)
  )
k_axis <- scale_x_continuous(breaks = c(2, seq(5, K_MAX, 5)), minor_breaks = NULL, limits = c(2, K_MAX),
                             expand = expansion(add = 0.6))
col_scale <- scale_colour_manual(values = COLS, labels = LABS, breaks = DIRS)
lt_scale <- scale_linetype_manual(values = c("Median" = "solid", "10th percentile" = "22"), drop = FALSE)
k_marks <- list(
  geom_vline(xintercept = BROAD_K, colour = "grey55", linewidth = 0.3, linetype = "dotted"),
  geom_vline(data = fine_k, aes(xintercept = K, colour = direction), linewidth = 0.3, linetype = "dotted",
             show.legend = FALSE)
)
mark_labels <- function(y) list(
  annotate("text", x = BROAD_K + 0.3, y = y, label = "broad", hjust = 0, vjust = 1, size = 2, colour = "grey40"),
  annotate("text", x = mean(fine_k$K) + 0.3, y = y, label = "fine", hjust = 0, vjust = 1, size = 2, colour = "grey40")
)

p_sil <- ggplot(sil_df, aes(K, mean_silhouette, colour = direction)) +
  k_marks +
  geom_line(linewidth = 0.45) + geom_point(size = 0.9) +
  mark_labels(max(sil_df$mean_silhouette) + 0.012) +
  k_axis + col_scale +
  scale_y_continuous(limits = c(0.3, max(sil_df$mean_silhouette) + 0.015), breaks = seq(0.3, 0.6, 0.05)) +
  labs(x = "Number of clusters (K)", y = "Mean silhouette width") +
  theme_pub

agreement_panel <- function(metric_name, ylab) {
  d <- filter(stab_long, metric == metric_name)
  ggplot(d, aes(K, value, colour = direction, linetype = stat)) +
    k_marks +
    geom_line(linewidth = 0.45) +
    geom_point(data = filter(d, stat == "Median"), size = 0.9, show.legend = FALSE) +
    mark_labels(1) +
    k_axis + col_scale + lt_scale +
    scale_y_continuous(limits = c(0.3, 1), breaks = seq(0.3, 1, 0.1)) +
    labs(x = "Number of clusters (K)", y = ylab) +
    guides(colour = "none") +  # direction legend comes from panel a
    theme_pub
}
p_ari <- agreement_panel("ARI", "Bootstrap ARI")
p_jac <- agreement_panel("Jaccard", "Mean cluster Jaccard")

cj <- bind_rows(lapply(DIRS, function(n) mutate(stab.ls[[n]]$cluster_summary, direction = n))) %>%
  filter(K %in% K_CLUSTER_PANEL) %>%
  mutate(cluster = factor(as.integer(cluster)),
         direction_lab = factor(LABS[direction], levels = LABS),
         K_lab = factor(paste0("K = ", K), levels = paste0("K = ", K_CLUSTER_PANEL)))
p_cl <- ggplot(cj, aes(cluster, median_Jaccard, colour = direction)) +
  geom_hline(yintercept = 0.5, colour = "grey80", linewidth = 0.25, linetype = "dashed") +
  geom_segment(aes(xend = cluster, y = q10_Jaccard, yend = median_Jaccard), linewidth = 0.5, show.legend = FALSE) +
  geom_point(size = 1.4, show.legend = FALSE) +
  facet_grid(direction_lab ~ K_lab, scales = "free_x", space = "free_x") +
  col_scale +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Cluster (full-data tree)", y = "Cluster Jaccard") +
  theme_pub +
  theme(panel.border = element_rect(fill = NA, colour = "grey70", linewidth = 0.3), axis.line = element_blank(),
        panel.spacing.x = unit(1.5, "mm"), strip.text.y = element_text(angle = 270))

top <- (p_sil | p_ari | p_jac)
fig <- (top / p_cl) +
  plot_layout(heights = c(1, 1.15), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1.5, "mm"), legend.spacing.x = unit(6, "mm"))

for (ext in c("pdf", "png")) {
  f <- file.path(FIGURE_DIR, paste0("clustering_metrics_composite.", ext))
  ggsave(f, fig, width = 180, height = 130, units = "mm", dpi = 300)
  message("wrote ", f)
}
cat("PLOT_COMPOSITE_OK\n")
