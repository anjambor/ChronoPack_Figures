suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]])))
} else {
  script_dir <- getwd()
}

# Copied from /new-stg/home/alex/tmp/Fig3/scripts/plot_dynamic_modes_broad_fine_examples_v3.R via Fig3_DynamicModes
# (its step 01, copied unchanged into v2 on 2026-11-16 so v2 builds both Fig3 figures; v1 is in ../old/).
# Inputs default to the ChronoPack r9 gene-trajectory tables in Fig3_DynamicModes/input/
# (INPUT_DIR / PC_DF_FILE / OUTPUT_DIR override, e.g. to rerun on the original Fig3 inputs).
project_dir <- dirname(script_dir)
input_dir <- Sys.getenv("INPUT_DIR", file.path(project_dir, "input"))
out_dir <- Sys.getenv("OUTPUT_DIR", file.path(project_dir, "output"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
args <- commandArgs(trailingOnly = TRUE)
cell_type_filter <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else NA_character_
output_suffix <- if (is.na(cell_type_filter)) "" else paste0("_", gsub("[^A-Za-z0-9]+", "_", cell_type_filter))
figure_version <- Sys.getenv("FIGURE_VERSION", "")
version_suffix <- if (nzchar(figure_version)) paste0("_", figure_version) else ""
fine_layout <- Sys.getenv("FINE_LAYOUT", "3x2")
if (!fine_layout %in% c("3x2", "2x2")) stop("FINE_LAYOUT must be '3x2' or '2x2'")
fine_ncol <- if (fine_layout == "2x2") 2 else 3
fine_slots <- 2 * fine_ncol
fine_grid_width <- 0.4 * fine_ncol
page_width <- 19.76 * (0.7 + fine_grid_width) / (0.7 + 1.2)
if (fine_layout == "2x2") version_suffix <- paste0(version_suffix, "_2x2")

pc_df_file <- Sys.getenv("PC_DF_FILE", "pc_df.new.csv")
message("Using pc_df: ", file.path(input_dir, pc_df_file))
pc_df <- read.csv(file.path(input_dir, pc_df_file), row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
traj_mat <- as.matrix(read.csv(file.path(input_dir, "gene_trajectory_fit.csv"), row.names = 1, check.names = FALSE))
tt <- read.csv(file.path(input_dir, "gene_trajectory_time_grid.csv"))$time

stopifnot(is.data.frame(pc_df))
stopifnot(is.matrix(traj_mat))
stopifnot(all(c(
  "gene",
  "cell_type",
  "broad_cluster",
  "cluster_fine_induced",
  "cluster_fine_repressed",
  "waldStat_1"
) %in% colnames(pc_df)))
stopifnot(length(tt) == ncol(traj_mat))

shared_rows <- intersect(rownames(traj_mat), rownames(pc_df))
if (!length(shared_rows)) stop("No shared rownames between traj_mat and pc_df.")

traj_mat <- traj_mat[shared_rows, , drop = FALSE]
pc_df <- pc_df[shared_rows, , drop = FALSE]

if (!is.na(cell_type_filter)) {
  keep_rows <- rownames(pc_df)[as.character(pc_df$cell_type) == cell_type_filter]
  if (!length(keep_rows)) {
    stop("No rows found for cell_type: ", cell_type_filter)
  }
  traj_mat <- traj_mat[keep_rows, , drop = FALSE]
  pc_df <- pc_df[keep_rows, , drop = FALSE]
}

z_row <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
}

X <- z_row(traj_mat)
colnames(X) <- paste0("t", seq_len(ncol(X)))

broad_order <- c(
  "impulse-up_1",
  "impulse-up_2",
  "trend-up",
  "trend-down_1",
  "trend-down_2",
  "transient-down"
)
broad_cols <- c(
  `impulse-up_1`   = "#D27286",
  `impulse-up_2`   = "#C25AA4",
  `trend-up`       = "#A5922A",
  `trend-down_1`   = "#4DA3C7",
  `trend-down_2`   = "#2F6FA3",
  `transient-down` = "#AF78DB"
)
induced_broad <- c("impulse-up_1", "impulse-up_2", "trend-up")
format_broad_label <- function(x) gsub("[-_]+", " ", x)

meta <- pc_df %>%
  mutate(
    row_id = rownames(pc_df),
    broad_cluster = as.character(broad_cluster),
    fine_cluster = if_else(
      broad_cluster %in% induced_broad,
      as.character(cluster_fine_induced),
      as.character(cluster_fine_repressed)
    ),
    fine_cluster_num = suppressWarnings(as.integer(fine_cluster)),
    broad_cluster = factor(broad_cluster, levels = broad_order),
    wald_rank = ifelse(is.finite(waldStat_1), waldStat_1, -Inf)
  ) %>%
  filter(!is.na(broad_cluster), !is.na(fine_cluster), fine_cluster != "-1")

if (!nrow(meta)) stop("No dynamic rows with valid v3 broad/fine cluster labels.")

time_lut <- tibble(tp = colnames(X), t = as.numeric(tt))

df_long <- as.data.frame(X, check.names = FALSE) %>%
  mutate(row_id = rownames(X)) %>%
  filter(row_id %in% meta$row_id) %>%
  pivot_longer(cols = -row_id, names_to = "tp", values_to = "z") %>%
  left_join(time_lut, by = "tp") %>%
  left_join(
    meta %>% select(row_id, gene, cell_type, fine_cluster, fine_cluster_num, broad_cluster),
    by = "row_id"
  )

summ_broad <- df_long %>%
  group_by(broad_cluster, t) %>%
  summarise(
    q25 = quantile(z, 0.25, na.rm = TRUE),
    q50 = quantile(z, 0.50, na.rm = TRUE),
    q75 = quantile(z, 0.75, na.rm = TRUE),
    n = n_distinct(row_id),
    .groups = "drop"
  )

summ_fine <- df_long %>%
  group_by(broad_cluster, fine_cluster, fine_cluster_num, t) %>%
  summarise(
    q25 = quantile(z, 0.25, na.rm = TRUE),
    q50 = quantile(z, 0.50, na.rm = TRUE),
    q75 = quantile(z, 0.75, na.rm = TRUE),
    n = n_distinct(row_id),
    .groups = "drop"
  )

example_df <- meta %>%
  arrange(broad_cluster, fine_cluster_num, desc(wald_rank), gene, cell_type) %>%
  group_by(broad_cluster, fine_cluster, fine_cluster_num) %>%
  distinct(gene, .keep_all = TRUE) %>%
  summarise(
    examples = paste(head(gene, 3), collapse = ", "),
    .groups = "drop"
  )

fine_label_df <- meta %>%
  count(broad_cluster, fine_cluster, fine_cluster_num, name = "n_gene_celltype") %>%
  arrange(broad_cluster, desc(n_gene_celltype), fine_cluster_num, fine_cluster) %>%
  left_join(example_df, by = c("broad_cluster", "fine_cluster", "fine_cluster_num")) %>%
  mutate(
    fine_cluster_label = paste0("Cluster ", fine_cluster, " (n=", n_gene_celltype, ")"),
    fine_label = paste0(
      fine_cluster_label,
      "\n",
      examples
    )
  )

write.csv(
  fine_label_df,
  file.path(out_dir, paste0("dynamic_modes_v3", output_suffix, "_fine_cluster_examples", version_suffix, ".csv")),
  row.names = FALSE
)

selected_df <- fine_label_df %>%
  group_by(broad_cluster) %>%
  arrange(desc(n_gene_celltype), fine_cluster_num, fine_cluster, .by_group = TRUE) %>%
  slice_head(n = 6) %>%
  ungroup()

write.csv(
  selected_df,
  file.path(out_dir, paste0("dynamic_modes_v3", output_suffix, "_selected_fine_cluster_examples", version_suffix, ".csv")),
  row.names = FALSE
)

plot_summary_curve <- function(df, title, color, subtitle = NULL, show_axes = TRUE, title_size = 11.2) {
  axis_theme <- if (show_axes) {
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 13)
    )
  } else {
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  }

  ggplot(df, aes(t, q50)) +
    geom_hline(yintercept = 0, linewidth = 0.22, color = "grey78") +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color, alpha = 0.20, color = NA) +
    geom_line(color = color, linewidth = 0.9) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), limits = range(tt)) +
    coord_cartesian(ylim = c(-2.1, 2.1)) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = title_size, lineheight = 0.95),
      plot.subtitle = element_text(size = 9.5, lineheight = 0.95),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    axis_theme +
    labs(
      title = title,
      subtitle = subtitle,
      x = if (show_axes) "Normalized Time" else NULL,
      y = if (show_axes) "z-score" else NULL
    )
}

row_plots <- lapply(broad_order, function(bc) {
  color <- broad_cols[[bc]]
  ids <- selected_df %>%
    filter(as.character(broad_cluster) == bc) %>%
    pull(fine_cluster)
  if (length(ids) > fine_slots) {
    stop(bc, " has ", length(ids), " fine clusters; FINE_LAYOUT=", fine_layout, " fits ", fine_slots)
  }
  # Four fine clusters in the 3x2 layout: lay out as 2 x 2 (as the original trend-up row did)
  ids_plot <- if (fine_ncol == 3 && length(ids) == 4) {
    c(ids[1:2], NA_character_, ids[3:4], NA_character_)
  } else {
    ids
  }

  broad_panel <- plot_summary_curve(
    df = summ_broad %>% filter(as.character(broad_cluster) == bc),
    title = format_broad_label(bc),
    color = color,
    title_size = 14.2
  )

  fine_panels <- lapply(ids_plot, function(fid) {
    if (is.na(fid)) {
      return(plot_spacer())
    }

    lab_row <- selected_df %>%
      filter(as.character(broad_cluster) == bc, fine_cluster == fid) %>%
      slice_head(n = 1)

    if (!nrow(lab_row) || is.na(lab_row$fine_cluster_label[[1]])) {
      return(plot_spacer())
    }

    plot_summary_curve(
      df = summ_fine %>% filter(as.character(broad_cluster) == bc, fine_cluster == fid),
      title = lab_row$fine_cluster_label[[1]],
      subtitle = lab_row$examples[[1]],
      color = color,
      show_axes = FALSE,
      title_size = 10.7
    )
  })

  blanks_needed <- fine_slots - length(fine_panels)
  if (blanks_needed > 0) {
    fine_panels <- c(fine_panels, rep(list(plot_spacer()), blanks_needed))
  }

  fine_grid <- wrap_plots(fine_panels, ncol = fine_ncol)
  broad_panel + fine_grid + plot_layout(widths = c(0.7, fine_grid_width))
})

single_page <- wrap_plots(row_plots[c(1, 4, 2, 5, 3, 6)], ncol = 2)

ggsave(
  file.path(out_dir, paste0("dynamic_modes_v3", output_suffix, "_broad_fine_examples", version_suffix, ".pdf")),
  single_page,
  width = page_width,
  height = 9.6
)
ggsave(
  file.path(out_dir, paste0("dynamic_modes_v3", output_suffix, "_broad_fine_examples", version_suffix, ".png")),
  single_page,
  width = page_width,
  height = 9.6,
  dpi = 300
)

message("Wrote:")
message("  ", file.path(out_dir, paste0("dynamic_modes_v3", output_suffix, "_broad_fine_examples", version_suffix, ".pdf")))
message("  ", file.path(out_dir, paste0("dynamic_modes_v3", output_suffix, "_broad_fine_examples", version_suffix, ".png")))
message("  ", file.path(out_dir, paste0("dynamic_modes_v3", output_suffix, "_fine_cluster_examples", version_suffix, ".csv")))
message("  ", file.path(out_dir, paste0("dynamic_modes_v3", output_suffix, "_selected_fine_cluster_examples", version_suffix, ".csv")))
