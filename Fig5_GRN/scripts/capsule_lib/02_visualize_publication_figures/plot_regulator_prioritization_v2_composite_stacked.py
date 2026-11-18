# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

from plot_regulator_prioritization_v2 import (
    FIG_DIR,
    METHOD_LABELS as STACKED_METHOD_LABELS,
    TF_COLORS,
    all_requested_tfs,
    requested_tfs,
)
from plot_regulator_prioritization_v2_nes_dot_comparison import (
    K562_METHODS,
    METHOD_LABELS as NES_METHOD_LABELS,
    SUMMARY_OUT,
    format_nes,
    load_k562_matrix,
    log_tick_formatter,
)


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables"
TOPK_PATH = TABLE_DIR / "regulator_prioritization_v2_fixed_topk_tf_proportions_full.csv"
OUT_STEM = "regulator_prioritization_v2_nes_composition_stacked_composite_signaling_only"
OUT_STEM_TF_ROWS = "regulator_prioritization_v2_nes_composition_stacked_composite_signaling_only_tf_rows"
# Fig5_GRN_v3 (2026-10-22, user request): the composite's two panels are saved as separate figures instead
OUT_STEM_LEFT = "regulator_prioritization_v2_nes_heatmap_signaling_only_tf_rows"  # rank-enrichment heatmap + fold
OUT_STEM_RIGHT = "regulator_prioritization_v2_top10000_composition_stacked_signaling_only"  # stacked bars
PANEL_B_CHRONO_METHODS = ["ChronoPack_K562", "ChronoPack_THP1"]
PANEL_B_COMPARATOR_METHODS = fig5_paths.insert_prior(["GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"])  # Fig5_GRN_v3
MILD_VIRIDIS = LinearSegmentedColormap.from_list(
    "mild_viridis",
    plt.get_cmap("viridis")(np.linspace(0.0, 0.90, 256)),
)


def text_color(log_value: float, norm: Normalize) -> str:
    r, g, b, _ = MILD_VIRIDIS(norm(log_value))
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "#1f1f1f" if luminance > 0.55 else "white"


def load_composition_pivot() -> pd.DataFrame:
    data = pd.read_csv(TOPK_PATH)
    sub = data[data["top_k"].eq(10000)].copy()
    order = all_requested_tfs()
    sub = sub[sub["is_evaluated_for_method"]].copy()
    plot = (
        sub.groupby(["method", "TF"], as_index=False)
        .agg(edge_share=("edge_share", "sum"))
    )
    comparator_totals = (
        plot[plot["method"].isin(PANEL_B_COMPARATOR_METHODS)]
        .groupby("method")["edge_share"]
        .sum()
        .sort_values(ascending=False)
    )
    method_order = PANEL_B_CHRONO_METHODS + comparator_totals.index.tolist()
    return (
        plot.pivot_table(index="method", columns="TF", values="edge_share", aggfunc="sum", fill_value=0)
        .reindex(index=method_order, columns=order, fill_value=0)
    )


def draw_transposed_matrix(ax, cax, matrix: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    transposed = matrix.T
    log_values = np.log10(transposed.to_numpy(float))
    finite = log_values[np.isfinite(log_values)]
    norm = Normalize(vmin=float(finite.min()), vmax=float(finite.max()))
    im = ax.imshow(log_values, cmap=MILD_VIRIDIS, norm=norm, aspect="auto")

    n_rows, n_cols = transposed.shape
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(transposed.columns.tolist(), rotation=0, fontsize=8.6)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels([NES_METHOD_LABELS[m] for m in transposed.index], fontsize=9.5)
    ax.tick_params(axis="both", length=0)
    ax.xaxis.tick_top()

    for i, method in enumerate(transposed.index):
        for j, tf in enumerate(transposed.columns):
            value = transposed.loc[method, tf]
            ax.text(
                j,
                i,
                format_nes(value),
                ha="center",
                va="center",
                fontsize=9.4,
                color=text_color(np.log10(value), norm),
                fontweight="bold" if method == "ChronoPack_K562" else "normal",
            )

    for x in np.arange(-0.5, n_cols, 1):
        ax.axvline(x, color="white", lw=1.0)
    for y in np.arange(-0.5, n_rows, 1):
        ax.axhline(y, color="white", lw=1.0)

    ax.axhline(0.5, color="#2b2b2b", lw=1.1)
    ax.add_patch(Rectangle((-0.5, -0.5), n_cols, 1, fill=False, edgecolor="#2b2b2b", linewidth=1.1))
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.spines[:].set_visible(False)

    fold_by_tf = summary.set_index("TF")["Chrono_vs_best_comparator_fold"]
    for j, tf in enumerate(transposed.columns):
        ax.text(
            j,
            n_rows + 0.10,
            f"{fold_by_tf.loc[tf]:.1f}$\\times$",
            ha="center",
            va="top",
            fontsize=7.8,
            fontweight="bold",
            color="#1f1f1f",
            clip_on=False,
        )
    ax.text(
        -0.65,
        n_rows + 0.10,
        "Chrono /\nbest comparator",
        ha="right",
        va="top",
        fontsize=7.8,
        color="#333333",
        clip_on=False,
    )

    cbar = plt.colorbar(im, cax=cax, orientation="horizontal")
    ticks = np.array([0.01, 0.1, 1, 10, 100])
    ticks = ticks[(np.log10(ticks) >= norm.vmin) & (np.log10(ticks) <= norm.vmax)]
    cbar.set_ticks(np.log10(ticks))
    cbar.set_ticklabels([log_tick_formatter(t, 0) for t in ticks])
    cbar.set_label("Rank-enrichment score (log scale)", fontsize=8.6, labelpad=2)
    cbar.ax.tick_params(labelsize=7.8, length=2, pad=1)
    return transposed


def draw_tf_rows_matrix(ax, cax, matrix: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    log_values = np.log10(matrix.to_numpy(float))
    finite = log_values[np.isfinite(log_values)]
    norm = Normalize(vmin=float(finite.min()), vmax=float(finite.max()))
    im = ax.imshow(log_values, cmap=MILD_VIRIDIS, norm=norm, aspect="auto")

    n_rows, n_cols = matrix.shape
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([NES_METHOD_LABELS[m] for m in matrix.columns], rotation=35, ha="right", fontsize=8.8)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(matrix.index.tolist(), fontsize=9.6)
    ax.tick_params(axis="both", length=0)
    ax.xaxis.tick_top()

    for i, tf in enumerate(matrix.index):
        for j, method in enumerate(matrix.columns):
            value = matrix.loc[tf, method]
            ax.text(
                j,
                i,
                format_nes(value),
                ha="center",
                va="center",
                fontsize=9.2,
                color=text_color(np.log10(value), norm),
                fontweight="bold" if method == "ChronoPack_K562" else "normal",
            )

    for x in np.arange(-0.5, n_cols, 1):
        ax.axvline(x, color="white", lw=1.0)
    for y in np.arange(-0.5, n_rows, 1):
        ax.axhline(y, color="white", lw=1.0)

    ax.axvline(0.5, color="#2b2b2b", lw=1.1)
    ax.add_patch(Rectangle((-0.5, -0.5), 1, n_rows, fill=False, edgecolor="#2b2b2b", linewidth=1.1))

    fold_by_tf = summary.set_index("TF")["Chrono_vs_best_comparator_fold"]
    fold_x = n_cols + 0.22
    for i, tf in enumerate(matrix.index):
        ax.text(
            fold_x,
            i,
            f"{fold_by_tf.loc[tf]:.1f}$\\times$",
            ha="left",
            va="center",
            fontsize=8.4,
            fontweight="bold",
            color="#1f1f1f",
            clip_on=False,
        )
    ax.text(
        fold_x,
        -0.95,
        "Chrono /\nbest comparator",
        ha="left",
        va="bottom",
        fontsize=8.2,
        color="#333333",
        clip_on=False,
    )

    ax.set_xlim(-0.5, n_cols + 1.35)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.spines[:].set_visible(False)

    cbar = plt.colorbar(im, cax=cax, orientation="vertical")
    ticks = np.array([0.01, 0.1, 1, 10, 100])
    ticks = ticks[(np.log10(ticks) >= norm.vmin) & (np.log10(ticks) <= norm.vmax)]
    cbar.set_ticks(np.log10(ticks))
    cbar.set_ticklabels([log_tick_formatter(t, 0) for t in ticks])
    cbar.set_label("Rank-enrichment score (log scale)", fontsize=8.6, labelpad=6)
    cbar.ax.tick_params(labelsize=7.8, length=2, pad=1)
    return matrix


def draw_row_aligned_stacked_bars(ax, pivot: pd.DataFrame) -> None:
    method_order = pivot.index.tolist()
    y = np.arange(len(method_order))
    left = np.zeros(len(method_order))
    order = all_requested_tfs()
    for tf in order:
        vals = pivot[tf].to_numpy(float)
        bars = ax.barh(
            y,
            vals,
            left=left,
            height=0.58,
            color=TF_COLORS.get(tf, "#BDBDBD"),
            edgecolor="white",
            linewidth=0.35,
        )
        for i, (bar, val) in enumerate(zip(bars, vals)):
            method = method_order[i]
            if tf in requested_tfs(method) and val >= 0.025:  # Fig5_GRN_v3: 0.018 let 'NFKB1' overflow a pySCENIC segment
                ax.text(
                    left[i] + val / 2,
                    y[i],
                    tf,
                    ha="center",
                    va="center",
                    fontsize=6.2,
                    color="white",
                )
        left += vals

    totals = pivot[order].sum(axis=1).to_numpy(float)
    x_max = max(0.04, float(totals.max()) * 1.18)
    label_pad = max(x_max * 0.018, 0.006)
    for yi, total in zip(y, totals):
        ax.text(
            total + label_pad,
            yi,
            f"{total * 100:.1f}%",
            ha="left",
            va="center",
            fontsize=8.2,
            fontweight="bold",
            color="#2a2a2a",
        )

    ax.set_yticks(y)
    ax.set_yticklabels([STACKED_METHOD_LABELS[m] for m in method_order], fontsize=9.0)
    ax.set_ylim(len(method_order) - 0.5, -0.5)
    ax.set_xlim(0, x_max)
    ax.set_xlabel("Share of top 10,000 predictions", fontsize=8.6, labelpad=8)
    ax.grid(axis="x", color="#e5e5e5", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


def draw_legend(fig, pivot: pd.DataFrame, x: float = 0.66, y: float = 0.04) -> None:
    legend_tfs = [tf for tf in all_requested_tfs() if pivot[tf].sum() > 0]
    handles = [plt.Rectangle((0, 0), 1, 1, color=TF_COLORS[tf]) for tf in legend_tfs]
    fig.legend(
        handles,
        legend_tfs,
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(x, y),
        ncol=6,
        fontsize=7.6,
        handlelength=1.2,
        columnspacing=0.9,
    )


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    matrix, summary = load_k562_matrix()
    summary.to_csv(SUMMARY_OUT, index=False)
    composition = load_composition_pivot()

    # Left panel: the composite's gridspec columns 0-1 (heatmap + colour bar), same panel size and margins
    fig = plt.figure(figsize=(13.0 * (0.98 + 0.035 + 0.06) / (0.98 + 0.035 + 0.24 + 1.05) + 1.4, 6.25))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.98, 0.035], left=0.16, right=0.86, top=0.94, bottom=0.18, wspace=0.06)
    draw_tf_rows_matrix(fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), matrix, summary)
    save_panel(fig, OUT_STEM_LEFT)

    # Right panel: the composite's gridspec column 3 (stacked bars) with the TF legend below it
    fig = plt.figure(figsize=(13.0 * 1.05 / (0.98 + 0.035 + 0.24 + 1.05) + 1.4, 6.25))
    ax_bars = fig.add_axes([0.2, 0.2, 0.76, 0.76])
    draw_row_aligned_stacked_bars(ax_bars, composition)
    draw_legend(fig, composition, x=0.58, y=0.035)
    save_panel(fig, OUT_STEM_RIGHT)


def save_panel(fig, stem: str) -> None:
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(FIG_DIR / f"{stem}.pdf")
    print(FIG_DIR / f"{stem}.png")


if __name__ == "__main__":
    main()
