# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns

import benchmark_grn_systems as b
from generate_system_broad_filtered_assets import BENCHMARK_LABELS
from generate_system_broad_filtered_publication_figures import VIOLIN_METHOD_LABELS
from generate_system_broad_filtered_publication_figures_v2 import ordered_tfs_for_benchmark
from plot_publication_v2_bounded_violins import format_outlier_value
from plot_publication_v2_violin_sepr import add_sepr


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables" / "publication_figures"
OUT_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "5b"

METRICS_TABLE = TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv"
METHOD_ORDER = b.METHOD_ORDER
A_PRIOR = 0.5
METRICS = [("sEPR", "sEPR"), ("AUPRC_ratio", "AUPRCr")]
BENCHMARK_XLIMS = {
    "TRRUST": {
        "sEPR": (0.0, 7.5),
        "AUPRC_ratio": (0.8, 8.0),
    },
    "K562_ChIP": {
        "sEPR": (0.9, 2.3),
        "AUPRC_ratio": (0.9, 2.0),
    },
}


def save(fig: plt.Figure, stem: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_metric(
    ax,
    data: pd.DataFrame,
    benchmark: str,
    metric: str,
    label: str,
    clipped_rows: list[dict],
    summary_rows: list[dict],
) -> None:
    method_labels = [VIOLIN_METHOD_LABELS[m] for m in METHOD_ORDER]
    tf_order = ordered_tfs_for_benchmark(benchmark)
    tf_offsets = {tf: offset for tf, offset in zip(tf_order, np.linspace(-0.32, 0.32, len(tf_order)))}
    xmin, configured_xmax = BENCHMARK_XLIMS[benchmark][metric]
    finite_max = float(data[metric].replace([np.inf, -np.inf], np.nan).max())
    xmax = configured_xmax if configured_xmax is not None else max(xmin + 0.1, finite_max * 1.25)
    x_range = xmax - xmin
    x_pad = x_range * 0.018
    outer_pad = x_range * 0.075

    plot_metric_col = f"{metric}_plot"
    data = data.copy()
    data[plot_metric_col] = data[metric].clip(lower=xmin, upper=configured_xmax)

    sns.violinplot(
        data=data,
        y="method_label",
        x=plot_metric_col,
        hue="method_label",
        order=method_labels,
        hue_order=method_labels,
        palette={VIOLIN_METHOD_LABELS[m]: b.METHOD_COLORS[m] for m in METHOD_ORDER},
        inner=None,
        cut=0,
        linewidth=1.0,
        saturation=0.65,
        density_norm="width",
        dodge=False,
        legend=False,
        ax=ax,
    )

    if benchmark == "TRRUST" and metric == "AUPRC_ratio":
        ax.set_xlim(xmin - outer_pad, xmax + outer_pad)
    else:
        ax.set_xlim(xmin, xmax)
    ax.axvline(xmin, color="#777777", lw=0.7, alpha=0.35)
    if configured_xmax is not None:
        ax.axvline(xmax, color="#777777", lw=0.7, alpha=0.35)
    ax.axvline(1, color="#333333", lw=0.8, ls="--", alpha=0.45)

    for method_i, method in enumerate(METHOD_ORDER):
        method_rows = data[data["method"] == method]
        for _, row in method_rows.iterrows():
            true_x = float(row[metric])
            plot_x = float(row[plot_metric_col])
            is_clipped = true_x < xmin or (configured_xmax is not None and true_x > xmax)
            y = method_i + tf_offsets[row["TF"]]
            ax.scatter(
                plot_x,
                y,
                s=23.4 if is_clipped else 19.5,
                color=b.METHOD_COLORS[method],
                edgecolor="black" if is_clipped else "white",
                linewidth=0.45,
                clip_on=False,
                zorder=4,
            )
            if is_clipped:
                clipped_rows.append(
                    {
                        "benchmark": benchmark,
                        "method": method,
                        "TF": row["TF"],
                        "metric": label,
                        "true_value": true_x,
                        "plotted_value": plot_x,
                    }
                )
                if true_x > xmax:
                    text_x = xmax + x_pad
                    ha = "left"
                else:
                    text_x = xmin - x_pad
                    ha = "right"
                ax.text(
                    text_x,
                    y,
                    f"{row['TF']} ({format_outlier_value(true_x)})",
                    ha=ha,
                    va="center",
                    fontsize=6.63,
                    color="#111111",
                    clip_on=False,
                    zorder=5,
                )
            else:
                ax.text(
                    plot_x + x_pad,
                    y,
                    row["TF"],
                    ha="left",
                    va="center",
                    fontsize=6.76,
                    color="#222222",
                    clip_on=False,
                    zorder=5,
                )

        vals = method_rows[metric].replace([np.inf, -np.inf], np.nan).dropna()
        if vals.empty:
            continue
        median_val = float(vals.median())
        median_plot = float(np.clip(median_val, xmin, configured_xmax))
        summary_rows.append(
            {
                "benchmark": benchmark,
                "method": method,
                "metric": label,
                "statistic": "median",
                "value": median_val,
                "plotted_value": median_plot,
                "n_TFs": int(vals.shape[0]),
            }
        )
        ax.scatter(
            median_plot,
            method_i,
            marker="D",
            s=42,
            facecolor="white",
            edgecolor="#111111",
            linewidth=0.9,
            clip_on=False,
            zorder=7,
        )

    ax.set_title(label, fontsize=12, pad=8)
    ax.set_xlabel(label)
    ax.set_ylabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_benchmark(data: pd.DataFrame, benchmark: str) -> None:
    bench = data[data["benchmark"] == benchmark].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.6 * max(1.0, len(METHOD_ORDER) / 8)), sharey=True)  # Fig5_GRN_v3
    summary_rows: list[dict] = []
    clipped_rows: list[dict] = []
    for ax, (metric, label) in zip(axes, METRICS):
        plot_metric(ax, bench, benchmark, metric, label, clipped_rows, summary_rows)

    method_labels = [VIOLIN_METHOD_LABELS[m] for m in METHOD_ORDER]
    axes[0].set_yticks(range(len(method_labels)))
    axes[0].set_yticklabels(method_labels, fontsize=9)
    axes[1].legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="D",
                linestyle="",
                markerfacecolor="white",
                markeredgecolor="#111111",
                markersize=5.5,
                label="Median",
            )
        ],
        frameon=False,
        loc="lower right",
        fontsize=8,
    )
    fig.suptitle(f"{BENCHMARK_LABELS[benchmark]} per-TF performance (sEPR a={A_PRIOR:g})", fontsize=13, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    stem = f"publication_violin_{benchmark}_per_tf_sEPR_AUPRCr_mean_median_all_methods_a0p5"
    save(fig, stem)

    clipped = pd.DataFrame(clipped_rows)
    summary = pd.DataFrame(summary_rows)
    clipped.to_csv(OUT_DIR / f"{stem}_clipped_points.csv", index=False)
    summary.to_csv(OUT_DIR / f"{stem}_summary.csv", index=False)
    print(OUT_DIR / f"{stem}.pdf")
    if not clipped.empty:
        print(clipped.to_string(index=False))


def main() -> None:
    data = add_sepr(pd.read_csv(METRICS_TABLE), a=A_PRIOR)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_benchmark(data, benchmark)


if __name__ == "__main__":
    main()
