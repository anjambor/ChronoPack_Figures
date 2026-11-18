# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

import benchmark_grn_systems as b
from generate_system_broad_filtered_publication_figures import VIOLIN_METHOD_LABELS
from generate_system_broad_filtered_publication_figures_v2 import ordered_tfs_for_benchmark
from plot_publication_v2_bounded_violins import format_outlier_value
from plot_publication_v2_violin_sepr import add_sepr


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables" / "publication_figures"
OUT_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "5b"

METRICS_TABLE = TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv"
A_PRIOR = 0.5
METHOD_ORDER = b.METHOD_ORDER
BENCHMARK_CONFIGS = {
    "TRRUST": {
        "label": "TRRUST",
        "metrics": [
            ("sEPR", "sEPR", 0.0, 7.5, 2.0),
            ("AUPRC_ratio", "AUPRCr", 0.5, 5.5, 1.4),
        ],
    },
    "K562_ChIP": {
        "label": "K562 ChIP",
        "metrics": [
            ("sEPR", "sEPR", 0.9, 2.3, 1.55),
            ("AUPRC_ratio", "AUPRCr", 0.9, 2.0, 1.45),
        ],
    },
}
LABEL_TFS = {"RELA", "NFKB1", "JUN", "FOS"}
OUTLIER_LABEL_TFS = {"REL", "ELK1"}
LABEL_ALL_POINTS = True
LABEL_Y_NUDGE = {
    "RELA": -0.060,
    "NFKB1": 0.018,
    "JUN": 0.080,
    "FOS": 0.000,
    "REL": 0.000,
    "ELK1": 0.000,
}


def method_label(method: str) -> str:
    return VIOLIN_METHOD_LABELS[method].split("\n")[0]


def softened_color(hex_color: str, mix: float = 0.18) -> tuple[float, float, float]:
    rgb = np.array([int(hex_color[i : i + 2], 16) for i in (1, 3, 5)], dtype=float) / 255.0
    return tuple(rgb * (1 - mix) + mix)


def tf_offsets(tf_order: list[str]) -> dict[str, float]:
    return {tf: offset for tf, offset in zip(tf_order, np.linspace(-0.37, 0.37, len(tf_order)))}


def plot_metric(
    ax,
    data: pd.DataFrame,
    benchmark: str,
    metric: str,
    title: str,
    xmin: float,
    xmax: float,
    label_threshold: float,
    y_positions: dict[str, float],
    clipped_rows: list[dict],
    summary_rows: list[dict],
) -> None:
    tf_order = ordered_tfs_for_benchmark(benchmark)
    offsets = tf_offsets(tf_order)
    x_range = xmax - xmin
    label_pad = x_range * 0.020
    edge_pad = x_range * 0.018

    ax.axvline(1, color="#9a9a9a", lw=0.9, ls="--", alpha=0.8, zorder=0)
    ax.grid(axis="x", color="#e7e7e7", linewidth=0.65)
    ax.set_axisbelow(True)

    for method in METHOD_ORDER:
        method_rows = data[data["method"].eq(method)].copy()
        if method_rows.empty:
            continue
        row_y = y_positions[method]
        values = method_rows[metric].replace([np.inf, -np.inf], np.nan).dropna()
        clipped_values = values.clip(lower=xmin, upper=xmax)
        if not values.empty:
            q1, q3 = np.percentile(values, [25, 75])
            ax.hlines(
                row_y,
                max(xmin, float(q1)),
                min(xmax, float(q3)),
                color="#2a2a2a",
                linewidth=1.35,
                alpha=0.22,
                zorder=1,
            )
            ax.hlines(
                row_y,
                float(clipped_values.min()),
                float(clipped_values.max()),
                color="#2a2a2a",
                linewidth=0.85,
                alpha=0.13,
                zorder=1,
            )

            median_val = float(values.median())
            median_plot = float(np.clip(median_val, xmin, xmax))
            summary_rows.append(
                {
                    "benchmark": benchmark,
                    "method": method,
                    "metric": title,
                    "statistic": "median",
                    "value": median_val,
                    "plotted_value": median_plot,
                    "n_TFs": int(values.shape[0]),
                }
            )
            ax.scatter(
                median_plot,
                row_y,
                marker="D",
                s=82,
                facecolor="white",
                edgecolor="#111111",
                linewidth=1.2,
                zorder=6,
            )

        color = softened_color(b.METHOD_COLORS[method], mix=0.10 if "ChronoPack" in method else 0.22)
        for _, row in method_rows.iterrows():
            true_x = float(row[metric])
            if not np.isfinite(true_x):
                continue
            tf = row["TF"]
            plot_x = float(np.clip(true_x, xmin, xmax))
            y = row_y + offsets.get(tf, 0.0)
            high_clip = true_x > xmax
            low_clip = true_x < xmin
            marker = ">" if high_clip else "<" if low_clip else "o"
            ax.scatter(
                plot_x,
                y,
                s=64 if "ChronoPack" in method else 50,
                marker=marker,
                color=color,
                edgecolor="#111111" if high_clip or low_clip else "white",
                linewidth=0.65,
                clip_on=False,
                zorder=4,
            )

            should_label = LABEL_ALL_POINTS or (
                (
                    high_clip
                    and tf in OUTLIER_LABEL_TFS
                ) or (
                    tf in LABEL_TFS
                    and ("ChronoPack" in method or true_x >= label_threshold)
                )
            )
            if not should_label:
                continue

            label_y = y + LABEL_Y_NUDGE.get(tf, 0.0)
            if high_clip:
                clipped_rows.append(
                    {
                        "benchmark": benchmark,
                        "method": method,
                        "TF": tf,
                        "metric": title,
                        "true_value": true_x,
                        "plotted_value": plot_x,
                    }
                )
                ax.annotate(
                    f"{tf} {format_outlier_value(true_x)}",
                    xy=(xmax, y),
                    xytext=(xmax + edge_pad, label_y),
                    ha="left",
                    va="center",
                    fontsize=8.2,
                    color="#222222",
                    arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.55},
                    clip_on=False,
                    zorder=8,
                )
            elif low_clip:
                clipped_rows.append(
                    {
                        "benchmark": benchmark,
                        "method": method,
                        "TF": tf,
                        "metric": title,
                        "true_value": true_x,
                        "plotted_value": plot_x,
                    }
                )
            else:
                ax.text(
                    plot_x + label_pad,
                    label_y,
                    tf,
                    ha="left",
                    va="center",
                    fontsize=7.6,
                    color="#2d2d2d",
                    clip_on=False,
                    zorder=5,
                )

    ax.set_title("")
    ax.set_xlabel(title, fontsize=12.2)
    ax.set_xlim(xmin, xmax)
    if metric in {"AUPRC_ratio", "EPR"} and xmin < 1.0:
        ax.set_xticks(np.arange(1.0, xmax + 0.001, 1.0))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=10.5)


def save(fig: plt.Figure, out_stem: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{out_stem}.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{out_stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_benchmark(all_data: pd.DataFrame, benchmark: str) -> None:
    config = BENCHMARK_CONFIGS[benchmark]
    data = all_data[all_data["benchmark"].eq(benchmark) & all_data["method"].isin(METHOD_ORDER)].copy()
    y_values, fig_height, (k562_label_y, thp1_label_y) = b.strip_rows(METHOD_ORDER, 4.65)  # Fig5_GRN_v3
    y_positions = dict(zip(METHOD_ORDER, y_values))
    y_ticks = [y_positions[m] for m in METHOD_ORDER]
    y_labels = [method_label(m) for m in METHOD_ORDER]

    fig, axes = plt.subplots(1, 2, figsize=(10.6, fig_height), sharey=True)
    clipped_rows: list[dict] = []
    summary_rows: list[dict] = []
    for ax, (metric, title, xmin, xmax, label_threshold) in zip(axes, config["metrics"]):
        plot_metric(
            ax,
            data,
            benchmark,
            metric,
            title,
            xmin,
            xmax,
            label_threshold,
            y_positions,
            clipped_rows,
            summary_rows,
        )

    for ax in axes:
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_labels, fontsize=11.0)
        ax.set_ylim(max(y_values) + 0.55, -0.62)
    axes[0].text(
        -0.30,
        k562_label_y,
        "K562",
        transform=axes[0].transAxes,
        ha="center",
        va="center",
        rotation=90,
        fontsize=9.8,
        color="#555555",
    )
    axes[0].text(
        -0.30,
        thp1_label_y,
        "THP1",
        transform=axes[0].transAxes,
        ha="center",
        va="center",
        rotation=90,
        fontsize=9.8,
        color="#555555",
    )

    axes[1].legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="D",
                linestyle="",
                markerfacecolor="white",
                markeredgecolor="#111111",
                markersize=7.8,
                label="Median",
            )
        ],
        frameon=False,
        loc="lower right",
        fontsize=10.0,
        handletextpad=0.4,
    )
    fig.suptitle(config["label"], fontsize=16.2, y=0.985)
    fig.tight_layout(rect=[0.055, 0.02, 1, 0.945], w_pad=2.0)

    out_stem = f"publication_strip_forest_{benchmark}_per_tf_sEPR_AUPRCr_all_methods_all_labeled_v2_a0p5"
    save(fig, out_stem)

    clipped = pd.DataFrame(clipped_rows)
    summary = pd.DataFrame(summary_rows)
    clipped.to_csv(OUT_DIR / f"{out_stem}_clipped_points.csv", index=False)
    summary.to_csv(OUT_DIR / f"{out_stem}_summary.csv", index=False)
    print(OUT_DIR / f"{out_stem}.pdf")
    print(OUT_DIR / f"{out_stem}.png")
    if not clipped.empty:
        print(clipped.to_string(index=False))


def main() -> None:
    data = add_sepr(pd.read_csv(METRICS_TABLE), a=A_PRIOR)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_benchmark(data, benchmark)


if __name__ == "__main__":
    main()
