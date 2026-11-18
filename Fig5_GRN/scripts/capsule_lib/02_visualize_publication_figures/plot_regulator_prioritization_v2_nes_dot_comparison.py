# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
import numpy as np
import pandas as pd


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables"
FIG_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "regulator_prioritization_v2"
METRICS_PATH = TABLE_DIR / "regulator_prioritization_v2_tf_rank_enrichment.csv"

OUT_STEM = "tf_rank_enrichment_nes_dot_comparison_v2"
SUMMARY_OUT = TABLE_DIR / "regulator_prioritization_v2_chrono_vs_best_comparator_nes.csv"

K562_METHODS = fig5_paths.insert_prior(["ChronoPack_K562", "GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"])  # Fig5_GRN_v3
COMPARATOR_METHODS = fig5_paths.insert_prior(["GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"])
K562_TFS = ["RELA", "NFKB1", "JUN", "ELK1", "CREB1", "ATF1", "FOS", "ATF2"]

METHOD_LABELS = {
    "ChronoPack_K562": "Chrono",
    "GENIE3": "GENIE3",
    "SINCERITIES": "SINCERITIES",
    "GRNBOOST2": "GRNBOOST2",
    "PPCOR": "PPCOR",
    **fig5_paths.PRIOR_LABELS,
}


def format_nes(value: float) -> str:
    if pd.isna(value):
        return ""
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.1f}"
    if value >= 1:
        return f"{value:.2f}".rstrip("0").rstrip(".")
    if value >= 0.01:
        return f"{value:.2f}"
    return f"{value:.3f}"


def log_tick_formatter(x: float, _pos: int) -> str:
    if x <= 0:
        return ""
    if x >= 100:
        return f"{x:.0f}"
    if x >= 1:
        return f"{x:g}"
    return f"{x:.3g}"


def load_k562_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    metrics = pd.read_csv(METRICS_PATH)
    metrics["NES"] = pd.to_numeric(metrics["NES"], errors="coerce")
    matrix = (
        metrics[metrics["method"].isin(K562_METHODS) & metrics["TF"].isin(K562_TFS)]
        .pivot_table(index="TF", columns="method", values="NES", aggfunc="first")
        .reindex(index=K562_TFS, columns=K562_METHODS)
    )
    summary = pd.DataFrame(
        {
            "TF": matrix.index,
            "Chrono_K562_NES": matrix["ChronoPack_K562"],
            "best_comparator_method": matrix[COMPARATOR_METHODS].idxmax(axis=1),
            "best_comparator_NES": matrix[COMPARATOR_METHODS].max(axis=1),
        }
    )
    summary["Chrono_vs_best_comparator_fold"] = summary["Chrono_K562_NES"] / summary["best_comparator_NES"]
    summary = summary.sort_values("Chrono_vs_best_comparator_fold", ascending=False).reset_index(drop=True)
    matrix = matrix.reindex(summary["TF"])
    return matrix, summary


def text_color(log_value: float, norm: Normalize) -> str:
    return "white" if norm(log_value) > 0.56 else "#1f1f1f"


def draw_group_headers(ax) -> None:
    transform = ax.get_xaxis_transform()
    y = 1.08
    ax.text(0, y, "OUR METHOD", transform=transform, ha="center", va="bottom", fontsize=8.5, color="#555555")
    ax.text((1 + len(COMPARATOR_METHODS)) / 2, y, "BASELINE METHODS", transform=transform, ha="center", va="bottom", fontsize=8.5, color="#555555")
    ax.plot([-0.45, 0.45], [1.065, 1.065], transform=transform, color="#777777", lw=0.8, clip_on=False)
    ax.plot([0.75, len(COMPARATOR_METHODS) + 0.25], [1.065, 1.065], transform=transform, color="#777777", lw=0.8, clip_on=False)


def plot_matrix(ax, cax, matrix: pd.DataFrame, summary: pd.DataFrame) -> None:
    log_values = np.log10(matrix.to_numpy(float))
    finite = log_values[np.isfinite(log_values)]
    norm = Normalize(vmin=float(finite.min()), vmax=float(finite.max()))
    im = ax.imshow(log_values, cmap="viridis", norm=norm, aspect="auto")

    n_rows, n_cols = matrix.shape
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([METHOD_LABELS[m] for m in matrix.columns], fontsize=8.4)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(matrix.index.tolist(), fontsize=10)
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
                fontsize=8.6,
                color=text_color(np.log10(value), norm),
                fontweight="bold" if method == "ChronoPack_K562" else "normal",
            )

    for x in np.arange(-0.5, n_cols, 1):
        ax.axvline(x, color="white", lw=1.0)
    for y in np.arange(-0.5, n_rows, 1):
        ax.axhline(y, color="white", lw=1.0)

    ax.axvline(0.5, color="#2b2b2b", lw=1.1)
    ax.add_patch(Rectangle((-0.5, -0.5), 1, n_rows, fill=False, edgecolor="#2b2b2b", linewidth=1.1))
    ax.set_xlim(-0.5, n_cols + 1.6)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_title("")
    ax.spines[:].set_visible(False)
    draw_group_headers(ax)

    fold_by_tf = summary.set_index("TF")["Chrono_vs_best_comparator_fold"]
    x_fold = n_cols + 0.9
    for i, tf in enumerate(matrix.index):
        ax.text(
            x_fold,
            i,
            f"{fold_by_tf.loc[tf]:.1f}$\\times$",
            ha="center",
            va="center",
            fontsize=9.4,
            fontweight="bold",
            color="#1f1f1f",
        )
    ax.text(
        x_fold,
        -0.74,
        "Chrono /\nbest comparator",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#333333",
    )
    ax.axvline(n_cols - 0.25, color="#d4d4d4", lw=0.9)

    cbar = plt.colorbar(im, cax=cax)
    ticks = np.array([0.01, 0.1, 1, 10, 100])
    ticks = ticks[(np.log10(ticks) >= norm.vmin) & (np.log10(ticks) <= norm.vmax)]
    cbar.set_ticks(np.log10(ticks))
    cbar.set_ticklabels([log_tick_formatter(t, 0) for t in ticks])
    cbar.ax.set_title("NES\n(log scale)", fontsize=8.2, pad=8)
    cbar.ax.tick_params(labelsize=8, length=2)


def plot_scatter(ax, summary: pd.DataFrame) -> None:
    x = summary["best_comparator_NES"].to_numpy(float)
    y = summary["Chrono_K562_NES"].to_numpy(float)
    lo = 1.5
    hi = 240
    ax.plot([lo, hi], [lo, hi], color="#9c9c9c", lw=0.9, ls="--", zorder=1)
    ax.scatter(x, y, s=42, color="#D99600", edgecolor="#111111", linewidth=0.6, zorder=3)
    label_positions = {
        "FOS": (3.2, 135),
        "ELK1": (6.9, 103),
        "RELA": (18.0, 207),
        "JUN": (12.8, 142),
        "CREB1": (33.0, 101),
        "ATF2": (32.0, 136),
        "NFKB1": (43.0, 175),
        "ATF1": (28.5, 62),
    }
    for _, row in summary.iterrows():
        text_xy = label_positions.get(row["TF"], (row["best_comparator_NES"], row["Chrono_K562_NES"]))
        ax.annotate(
            row["TF"],
            (row["best_comparator_NES"], row["Chrono_K562_NES"]),
            xytext=text_xy,
            textcoords="data",
            fontsize=7.9,
            color="#222222",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.5},
            arrowprops={"arrowstyle": "-", "color": "#9d9d9d", "lw": 0.45, "shrinkA": 2, "shrinkB": 2},
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.xaxis.set_major_formatter(FuncFormatter(log_tick_formatter))
    ax.yaxis.set_major_formatter(FuncFormatter(log_tick_formatter))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(axis="both", which="major", color="#dedede", lw=0.7)
    ax.grid(axis="both", which="minor", visible=False)
    ax.set_xlabel("Best comparator NES", fontsize=9)
    ax.set_ylabel("Chrono NES", fontsize=9)
    ax.set_title("Chrono vs best baseline", loc="left", fontsize=11, fontweight="bold", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    matrix, summary = load_k562_matrix()
    summary.to_csv(SUMMARY_OUT, index=False)
    print(summary.to_string(index=False, formatters={"Chrono_vs_best_comparator_fold": "{:.2f}".format}))

    fig = plt.figure(figsize=(12.0, 5.45))
    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[2.0, 0.055, 0.92],
        wspace=0.45,
        left=0.07,
        right=0.98,
        top=0.78,
        bottom=0.13,
    )
    ax_matrix = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    ax_scatter = fig.add_subplot(gs[0, 2])

    plot_matrix(ax_matrix, cax, matrix, summary)
    plot_scatter(ax_scatter, summary)
    fig.text(0.07, 0.94, "TF-target enrichment in K562", ha="left", va="top", fontsize=14, fontweight="bold")

    fig.savefig(FIG_DIR / f"{OUT_STEM}.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{OUT_STEM}.pdf", bbox_inches="tight")
    plt.close(fig)

    print(FIG_DIR / f"{OUT_STEM}.pdf")
    print(FIG_DIR / f"{OUT_STEM}.png")


if __name__ == "__main__":
    main()
