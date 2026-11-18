# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import benchmark_grn_systems as b
from generate_system_broad_filtered_assets import BENCHMARK_LABELS
from generate_system_broad_filtered_publication_figures import VIOLIN_METHOD_LABELS
from generate_system_broad_filtered_publication_figures_v2 import ALL_METHOD_ORDER, ordered_tfs_for_benchmark


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables" / "publication_figures"
OUT_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "5b"
PER_TF_TABLE = TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics.csv"

LIMITS = {
    "K562_ChIP": {
        "EPR": (1.0, 2.3),
        "AUPRC_ratio": (1.0, 2.0),
    },
    "TRRUST": {
        "EPR": (1.0, 15.0),
        "AUPRC_ratio": (1.0, 5.0),
    },
}

METRIC_LABELS = {
    "EPR": "EPR",
    "AUPRC_ratio": "AUPRCr",
}


def format_outlier_value(value: float) -> str:
    if abs(value) >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}"


def save(fig: plt.Figure, stem: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_bounded_violin(data: pd.DataFrame, benchmark: str) -> pd.DataFrame:
    data = data[data["benchmark"] == benchmark].copy()
    method_labels = [VIOLIN_METHOD_LABELS[m] for m in ALL_METHOD_ORDER]
    tf_order = ordered_tfs_for_benchmark(benchmark)
    tf_offsets = {tf: offset for tf, offset in zip(tf_order, np.linspace(-0.32, 0.32, len(tf_order)))}
    palette = {VIOLIN_METHOD_LABELS[m]: b.METHOD_COLORS[m] for m in ALL_METHOD_ORDER}

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.8), sharey=True)
    clipped_rows = []
    for ax, metric in zip(axes, ["EPR", "AUPRC_ratio"]):
        xmin, xmax = LIMITS[benchmark][metric]
        label = METRIC_LABELS[metric]
        plot_data = data.copy()
        plot_data[f"{metric}_plot"] = plot_data[metric].clip(lower=xmin, upper=xmax)
        sns.violinplot(
            data=plot_data,
            y="method_label",
            x=f"{metric}_plot",
            hue="method_label",
            order=method_labels,
            hue_order=method_labels,
            palette=palette,
            inner=None,
            cut=0,
            linewidth=1.0,
            saturation=0.65,
            density_norm="width",
            dodge=False,
            legend=False,
            ax=ax,
        )
        x_range = xmax - xmin
        x_pad = x_range * 0.018
        outer_pad = x_range * 0.075
        ax.set_xlim(xmin - outer_pad, xmax + outer_pad)
        ax.axvline(xmin, color="#777777", lw=0.7, alpha=0.35)
        ax.axvline(xmax, color="#777777", lw=0.7, alpha=0.35)
        ax.axvline(1, color="#333333", lw=0.8, ls="--", alpha=0.45)

        for method_i, method in enumerate(ALL_METHOD_ORDER):
            method_rows = plot_data[plot_data["method"] == method]
            for _, row in method_rows.iterrows():
                true_x = float(row[metric])
                plot_x = float(row[f"{metric}_plot"])
                is_clipped = true_x < xmin or true_x > xmax
                y = method_i + tf_offsets[row["TF"]]
                ax.scatter(
                    plot_x,
                    y,
                    s=18 if is_clipped else 16,
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
                    text = f"{row['TF']} {format_outlier_value(true_x)}"
                    ax.text(text_x, y, text, ha=ha, va="center", fontsize=5.3, color="#111111", clip_on=False, zorder=5)
                else:
                    ax.text(
                        plot_x + x_pad,
                        y,
                        row["TF"],
                        ha="left",
                        va="center",
                        fontsize=5.4,
                        color="#222222",
                        clip_on=False,
                        zorder=5,
                    )

        ax.set_title(label, fontsize=12, pad=8)
        ax.set_xlabel(label)
        ax.set_ylabel("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_yticks(range(len(method_labels)))
    axes[0].set_yticklabels(method_labels, fontsize=9)
    fig.suptitle(f"{BENCHMARK_LABELS[benchmark]} bounded axes", fontsize=13, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save(fig, f"publication_violin_{benchmark}_per_tf_EPR_AUPRCr_bounded")
    return pd.DataFrame(clipped_rows)


def main() -> None:
    data = pd.read_csv(PER_TF_TABLE)
    data = data[data["scope"] == "per_tf"].copy()
    data["EPR"] = data["EPR_raw"]
    data = data[["benchmark", "method", "method_label", "TF", "EPR", "AUPRC_ratio"]]
    clipped = []
    for benchmark in ["TRRUST", "K562_ChIP"]:
        clipped.append(plot_bounded_violin(data, benchmark))
    clipped_out = pd.concat(clipped, ignore_index=True)
    clipped_out.to_csv(TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_violin_bounded_clipped_points.csv", index=False)
    print(clipped_out.to_string(index=False))


if __name__ == "__main__":
    main()
