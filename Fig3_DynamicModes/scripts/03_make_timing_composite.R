suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# Fig3_DynamicModes_v2: copy of Fig3_DynamicModes/scripts/03_make_timing_composite.R (itself from
# /new-stg/home/alex/tmp/Fig3/scripts/make_publication_dynamic_modes_v10_vector.R). Mode selection,
# bars, panels and styling are unchanged; the downsampling fits are now ChronoPack fits on r9 data:
# curves, 95% bands and observed means come from output/chronopack/{curves,observed_means}.csv and
# onset/peak (true minutes) from output/chronopack/trajectory_descriptors.csv (steps 01-02).
script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]))) else getwd()
project_dir <- dirname(script_dir)

context <- list(
  context = "K562a_TNFa_10min_vs_30min",
  dir = file.path(project_dir, "output", "chronopack"),
  ref = "K562_trunc",
  ds = "K562_ds",
  cell_type = "K562a"
)

mode_df_path <- Sys.getenv("MODE_DF_PATH", file.path(project_dir, "input", "pc_df.new.csv"))
message("Using mode table: ", mode_df_path)
desc_path <- Sys.getenv("DESC_PATH", file.path(context$dir, "trajectory_descriptors.csv"))
message("Using timing descriptors: ", desc_path)
read_desc <- function(path) {
  if (grepl("\\.rds$", path)) readRDS(path) else read.csv(path, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
}
read_mode_df <- function(path) read.csv(path, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
mode_col <- "broad_cluster"
require_qvalue <- Sys.getenv("REQUIRE_QVALUE", unset = "false") == "true"
suffix_extra <- Sys.getenv("OUTPUT_SUFFIX_EXTRA", unset = "")
use_break_axis <- Sys.getenv("USE_BREAK_AXIS", unset = "true") != "false"
plot_modes <- c("impulse-up_1+2", "trend-up", "trend-down_1+2", "transient-down")
publication_genes <- c("NFKBIA", "CXCL1", "ATF3", "TNFAIP2", "IRF1", "HES1", "BIRC3", "CD69")
out_dir <- Sys.getenv("OUTPUT_DIR", file.path(project_dir, "output"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_base <- file.path(
  out_dir,
  paste0(
    "publication_timing_gene_composite_dynamic_modes_v10_original_dynamic_modes_v3_p005_fcgt05_all_assigned_combined_modes_with_transient",
    suffix_extra,
    "_vector_v2"
  )
)

break_axis_trans <- function(break_low = 1, break_high = 3) {
  gap <- break_high - break_low
  scales::trans_new(
    name = paste0("axis_break_", break_low, "_", break_high),
    transform = function(y) ifelse(y > break_low, y - gap, y),
    inverse = function(y) ifelse(y > break_low, y + gap, y)
  )
}

prepare_compare_df <- function() {
  desc <- read_desc(desc_path)
  mode_pc_df <- read_mode_df(mode_df_path)

  selected <- mode_pc_df %>%
    filter(
      cell_type == context$cell_type,
      !is.na(.data[[mode_col]]),
      .data[[mode_col]] != "low_significance",
      abs(peak_log2fc) > 0.5,
      pvalue_1 < 0.05,
      if (require_qvalue) qval_1 < 0.05 else TRUE
    ) %>%
    mutate(
      dynamic_mode = case_when(
        .data[[mode_col]] %in% c("impulse-up_1", "impulse-up_2") ~ "impulse-up_1+2",
        .data[[mode_col]] %in% c("trend-down_1", "trend-down_2") ~ "trend-down_1+2",
        TRUE ~ as.character(.data[[mode_col]])
      )
    ) %>%
    filter(dynamic_mode %in% plot_modes) %>%
    arrange(desc(waldStat_1)) %>%
    select(gene, dynamic_mode, waldStat_1, peak_log2fc)

  timing <- desc %>%
    filter(gene %in% selected$gene, cell_type %in% c(context$ref, context$ds)) %>%
    select(gene, cell_type, t_onset, t_peak) %>%
    pivot_wider(names_from = cell_type, values_from = c(t_onset, t_peak), names_sep = "__")

  selected %>%
    inner_join(timing, by = "gene") %>%
    mutate(
      onset_diff = .data[[paste0("t_onset__", context$ds)]] - .data[[paste0("t_onset__", context$ref)]],
      peak_diff = .data[[paste0("t_peak__", context$ds)]] - .data[[paste0("t_peak__", context$ref)]]
    )
}

check_strict_subset_of_original <- function(compare_df) {
  if (!require_qvalue) return(invisible(NULL))

  mode_pc_df <- read_mode_df(mode_df_path)
  original_genes <- mode_pc_df %>%
    filter(
      cell_type == context$cell_type,
      !is.na(.data[[mode_col]]),
      .data[[mode_col]] != "low_significance",
      abs(peak_log2fc) > 0.5,
      pvalue_1 < 0.05
    ) %>%
    mutate(
      dynamic_mode = case_when(
        .data[[mode_col]] %in% c("impulse-up_1", "impulse-up_2") ~ "impulse-up_1+2",
        .data[[mode_col]] %in% c("trend-down_1", "trend-down_2") ~ "trend-down_1+2",
        TRUE ~ as.character(.data[[mode_col]])
      )
    ) %>%
    filter(dynamic_mode %in% plot_modes) %>%
    pull(gene) %>%
    unique()

  new_genes <- unique(compare_df$gene)
  if (!all(new_genes %in% original_genes)) {
    stop("q-value filtered genes are not a subset of the original p-value filtered gene set")
  }
  if (length(new_genes) >= length(original_genes)) {
    stop("q-value filtered gene set is not a strict subset of the original p-value filtered gene set")
  }
  message("Subset check passed: ", length(new_genes), " q-filtered genes vs ", length(original_genes), " original p-filtered genes")
}

make_mode_bar <- function(df, value_col, title) {
  mode_labels <- c(
    `impulse-up_1+2` = "Impulse up 1+2",
    `trend-up` = "Trend up",
    `trend-down_1+2` = "Trend down 1+2",
    `transient-down` = "Transient down"
  )
  mode_colors <- c(
    `impulse-up_1+2` = "#CC668C",
    `trend-up` = "#A5922A",
    `trend-down_1+2` = "#3F88B5",
    `transient-down` = "#AF78DB"
  )
  y_break_low <- 1
  y_break_high <- 3

  summary_df <- df %>%
    filter(dynamic_mode %in% plot_modes, is.finite(.data[[value_col]])) %>%
    mutate(
      dynamic_mode = factor(dynamic_mode, levels = plot_modes),
      plot_value = abs(.data[[value_col]])
    ) %>%
    group_by(dynamic_mode) %>%
    summarise(
      n = n(),
      mean_value = mean(plot_value, na.rm = TRUE),
      se_value = sd(plot_value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      se_value = if_else(is.na(se_value), 0, se_value),
      errorbar_ymin = pmax(0, mean_value - se_value),
      errorbar_ymax = mean_value + se_value,
      label = paste0(sprintf("%.1f", mean_value), " +/- ", sprintf("%.1f", se_value))
    )
  if (use_break_axis) {
    summary_df <- summary_df %>%
      mutate(errorbar_ymin = if_else(errorbar_ymin > y_break_low & errorbar_ymin < y_break_high, y_break_high, errorbar_ymin))
  }

  y_max <- max(summary_df$mean_value + summary_df$se_value, na.rm = TRUE)
  y_lims <- c(0, y_max * 1.18)
  y_breaks <- c(0, y_break_low, seq(4, ceiling(y_lims[2] / 5) * 5, by = 5))

  p <- ggplot(summary_df, aes(x = dynamic_mode, y = mean_value, fill = dynamic_mode)) +
    geom_col(width = 0.68, color = "grey15", linewidth = 0.45, alpha = 0.82) +
    geom_errorbar(aes(ymin = errorbar_ymin, ymax = errorbar_ymax), width = 0.18, linewidth = 0.7, color = "black") +
    geom_text(aes(y = mean_value + se_value, label = label), nudge_y = 0.35, vjust = 0, size = 3.2, color = "black") +
    scale_x_discrete(labels = mode_labels, drop = FALSE, expand = expansion(add = 0.8)) +
    scale_fill_manual(values = mode_colors, guide = "none", drop = FALSE) +
    coord_cartesian(ylim = y_lims, clip = "off") +
    labs(title = title, x = NULL, y = "Absolute offset (min)") +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5),
      axis.title = element_text(size = 12.5),
      axis.text = element_text(color = "black", size = 11),
      axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1)
    )
  if (use_break_axis) {
    p <- p +
      annotate("text", x = 0.21, y = y_break_high + 0.08, label = "//", angle = 55, size = 4.3, color = "black") +
      scale_y_continuous(
        trans = break_axis_trans(y_break_low, y_break_high),
        breaks = y_breaks,
        expand = expansion(mult = c(0, 0.04))
      )
  } else {
    p <- p +
      scale_y_continuous(
        breaks = scales::pretty_breaks(n = 6),
        expand = expansion(mult = c(0, 0.04))
      )
  }
  p
}

role_of <- function(run) if_else(run == context$ref, "reference", "downsampled")

read_curves <- function() read.csv(file.path(context$dir, "curves.csv"), stringsAsFactors = FALSE)

make_curve_tbl <- function(genes) {
  read_curves() %>%
    filter(gene %in% genes, run %in% c(context$ref, context$ds)) %>%
    transmute(gene, cell_type = role_of(run), time_index = t, expression = mu_log1p)
}

make_observed_means_tbl <- function(genes) {
  read.csv(file.path(context$dir, "observed_means.csv"), stringsAsFactors = FALSE) %>%
    filter(gene %in% genes, run %in% c(context$ref, context$ds)) %>%
    transmute(gene, cell_type = role_of(run), observed_x = t, observed_expression = mean_log1p)
}

make_fit_ribbon_tbl <- function(genes) {
  read_curves() %>%
    filter(gene %in% genes, run %in% c(context$ref, context$ds)) %>%
    transmute(gene, cell_type = role_of(run), time_index = t, fit_lo = lo_log1p, fit_hi = hi_log1p)
}

make_timing_markers <- function(compare_df, genes) {
  ref <- context$ref
  ds <- context$ds
  bind_rows(
    compare_df %>% filter(gene %in% genes) %>% transmute(gene, cell_type = "reference", event = "Onset", event_time = .data[[paste0("t_onset__", ref)]]),
    compare_df %>% filter(gene %in% genes) %>% transmute(gene, cell_type = "reference", event = "Peak", event_time = .data[[paste0("t_peak__", ref)]]),
    compare_df %>% filter(gene %in% genes) %>% transmute(gene, cell_type = "downsampled", event = "Onset", event_time = .data[[paste0("t_onset__", ds)]]),
    compare_df %>% filter(gene %in% genes) %>% transmute(gene, cell_type = "downsampled", event = "Peak", event_time = .data[[paste0("t_peak__", ds)]])
  ) %>%
    filter(is.finite(event_time))
}

make_gene_panel <- function(curves, observed_means, ribbons, timing_markers, gene, show_y = FALSE, show_legend = FALSE) {
  plot_df <- curves %>% filter(gene == !!gene)
  obs_means <- observed_means %>% filter(gene == !!gene)
  ribbon_df <- ribbons %>% filter(gene == !!gene)
  marker_df <- bind_rows(lapply(seq_len(nrow(timing_markers %>% filter(gene == !!gene))), function(i) {
    marker <- (timing_markers %>% filter(gene == !!gene))[i, , drop = FALSE]
    curve_role <- plot_df %>% filter(cell_type == marker$cell_type)
    if (nrow(curve_role) == 0 || !is.finite(marker$event_time)) return(NULL)
    nearest <- which.min(abs(curve_role$time_index - marker$event_time))
    marker %>% mutate(marker_x = curve_role$time_index[nearest], marker_y = curve_role$expression[nearest])
  }))

  ggplot(plot_df, aes(x = time_index, y = expression, color = cell_type)) +
    geom_ribbon(
      data = ribbon_df,
      aes(x = time_index, ymin = fit_lo, ymax = fit_hi, fill = cell_type),
      inherit.aes = FALSE,
      alpha = 0.18,
      color = NA
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(data = obs_means, aes(x = observed_x, y = observed_expression, color = cell_type), inherit.aes = FALSE, shape = 21, fill = "white", size = 2.25, stroke = 0.65) +
    geom_point(data = marker_df, aes(x = marker_x, y = marker_y, shape = event), inherit.aes = FALSE, color = "firebrick", fill = "white", size = 3.1, stroke = 1.0) +
    scale_color_manual(values = c("reference" = "#386cb0", "downsampled" = "#222021"), labels = c("reference" = "Reference", "downsampled" = "Downsampled")) +
    scale_shape_manual(values = c("Onset" = 23, "Peak" = 24), breaks = c("Onset", "Peak")) +
    scale_fill_manual(values = c("reference" = "#386cb0", "downsampled" = "#222021"), guide = "none") +
    labs(title = gene, x = "Time (min)", y = if (show_y) "Expression" else NULL, color = NULL, shape = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text = element_text(color = "black", size = 10.5),
      axis.title = element_text(size = 12),
      legend.position = if (show_legend) "bottom" else "none",
      legend.text = element_text(size = 13),
      legend.key.width = unit(1.4, "cm"),
      legend.key.height = unit(0.5, "cm")
    )
}

message("Preparing v10 timing/mode summary")
compare_df <- prepare_compare_df()
check_strict_subset_of_original(compare_df)
write.csv(compare_df, paste0(out_base, "_genes.csv"), row.names = FALSE)
mode_summary <- bind_rows(lapply(c(peak = "peak_diff", onset = "onset_diff"), function(col) {
  compare_df %>%
    filter(dynamic_mode %in% plot_modes, is.finite(.data[[col]])) %>%
    mutate(dynamic_mode = factor(dynamic_mode, levels = plot_modes), v = abs(.data[[col]])) %>%
    group_by(dynamic_mode) %>%
    summarise(n = n(), mean_abs_offset_min = mean(v), se = sd(v) / sqrt(n()), .groups = "drop")
}), .id = "event")
write.csv(mode_summary, paste0(out_base, "_mode_summary.csv"), row.names = FALSE)
print(as.data.frame(mode_summary))
message("Loading trajectory curves")
curves <- make_curve_tbl(publication_genes)
message("Loading observed means")
observed_means <- make_observed_means_tbl(publication_genes)
message("Loading fit ribbons")
ribbons <- make_fit_ribbon_tbl(publication_genes)
message("Preparing timing markers")
timing_markers <- make_timing_markers(compare_df, publication_genes)

p_peak <- make_mode_bar(compare_df, "peak_diff", "Peak time absolute offsets by mode")
p_onset <- make_mode_bar(compare_df, "onset_diff", "Onset time absolute offsets by mode")
gene_panels <- lapply(seq_along(publication_genes), function(i) {
  make_gene_panel(
    curves,
    observed_means,
    ribbons,
    timing_markers,
    publication_genes[[i]],
    show_y = i %in% c(1, 5),
    show_legend = i == 8
  )
})

composite <- (
  p_peak + gene_panels[[1]] + gene_panels[[2]] + gene_panels[[3]] + gene_panels[[4]] +
    plot_layout(widths = c(1.87, 1, 1, 1, 1))
) / (
  p_onset + gene_panels[[5]] + gene_panels[[6]] + gene_panels[[7]] + gene_panels[[8]] +
    plot_layout(widths = c(1.87, 1, 1, 1, 1))
) / patchwork::guide_area() +
  plot_layout(heights = c(1, 1, 0.13), guides = "collect") +
  plot_annotation(
    title = "K562a TNF-alpha downsampling agreement",
    subtitle = paste0("Mean +/- SE of absolute timing offsets by ", mode_col, " dynamic response mode with representative fitted NB-GAM trajectories"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 19, hjust = 0),
      plot.subtitle = element_text(size = 14, hjust = 0),
      legend.position = "bottom",
      legend.text = element_text(size = 15),
      legend.key.width = unit(1.8, "cm"),
      legend.key.height = unit(0.65, "cm")
    )
  )

message("Writing PDF")
ggsave(paste0(out_base, ".pdf"), composite, device = cairo_pdf, width = 18.5, height = 7.2, units = "in", limitsize = FALSE)
message("Writing PNG preview")
ggsave(paste0(out_base, ".png"), composite, width = 18.5, height = 7.2, units = "in", dpi = 300, limitsize = FALSE)

message("Wrote ", paste0(out_base, ".pdf"))
message("Wrote ", paste0(out_base, ".png"))
