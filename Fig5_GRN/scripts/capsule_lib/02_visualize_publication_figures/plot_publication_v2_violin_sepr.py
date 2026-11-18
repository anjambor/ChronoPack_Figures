# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

import argparse
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
from plot_publication_v2_bounded_violins import LIMITS, format_outlier_value


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables" / "publication_figures"
OUT_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "5b"

METRICS_TABLE = TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv"
SEPR_TABLE = TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_sEPR_metrics_all_methods.csv"
UNBOUNDED_SUMMARY = TABLE_DIR / "system_broad_filtered_publication_v2_violin_sEPR_mean_median_all_methods.csv"
BOUNDED_SUMMARY = TABLE_DIR / "system_broad_filtered_publication_v2_bounded_violin_sEPR_mean_median_all_methods.csv"
BOUNDED_CLIPPED = TABLE_DIR / "system_broad_filtered_publication_v2_bounded_violin_sEPR_clipped_points_all_methods.csv"

METHOD_ORDER = b.METHOD_ORDER
METRICS = [("sEPR", "sEPR"), ("AUPRC_ratio", "AUPRCr")]
A_PRIOR = 1.0


def a_suffix(a: float) -> str:
    if np.isclose(a, A_PRIOR):
        return ""
    return f"_a{str(a).replace('.', 'p')}"


def suffixed_path(path: Path, a: float) -> Path:
    return path.with_name(f"{path.stem}{a_suffix(a)}{path.suffix}")


def save(fig: plt.Figure, stem: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def add_sepr(data: pd.DataFrame, a: float = A_PRIOR) -> pd.DataFrame:
    out = data[data["scope"] == "per_tf"].copy()
    out["K"] = pd.to_numeric(out["k_pos"], errors="coerce")
    out["M"] = pd.to_numeric(out["N"], errors="coerce")
    out["k"] = pd.to_numeric(out["k_EPR"], errors="coerce")
    out["X"] = pd.to_numeric(out["TP_at_k"], errors="coerce")
    out["p0"] = out["K"] / out["M"]
    out["lambda"] = out["k"] * out["p0"]
    out["EPR_raw"] = pd.to_numeric(out["EPR_raw"], errors="coerce")
    out["sEPR"] = (out["X"] + a) / (out["lambda"] + a)
    out["sEPR_a"] = a
    return out


def metric_limits(benchmark: str, metric: str) -> tuple[float, float]:
    if metric == "sEPR":
        return LIMITS[benchmark]["EPR"]
    return LIMITS[benchmark][metric]


def plot_points_and_stats(
    ax,
    data: pd.DataFrame,
    benchmark: str,
    metric: str,
    label: str,
    bounded: bool,
    clipped_rows: list[dict],
    summary_rows: list[dict],
) -> None:
    method_labels = [VIOLIN_METHOD_LABELS[m] for m in METHOD_ORDER]
    tf_order = ordered_tfs_for_benchmark(benchmark)
    tf_offsets = {tf: offset for tf, offset in zip(tf_order, np.linspace(-0.32, 0.32, len(tf_order)))}

    if bounded:
        xmin, xmax = metric_limits(benchmark, metric)
        plot_metric = f"{metric}_plot"
        data = data.copy()
        data[plot_metric] = data[metric].clip(lower=xmin, upper=xmax)
        x_range = xmax - xmin
        x_pad = x_range * 0.018
        outer_pad = x_range * 0.075
    else:
        plot_metric = metric
        xmin = 0.0
        finite_max = float(data[metric].replace([np.inf, -np.inf], np.nan).max())
        xmax = max(1.2, finite_max * 1.25)
        x_pad = xmax * 0.012

    sns.violinplot(
        data=data,
        y="method_label",
        x=plot_metric,
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

    if bounded:
        ax.set_xlim(xmin - outer_pad, xmax + outer_pad)
        ax.axvline(xmin, color="#777777", lw=0.7, alpha=0.35)
        ax.axvline(xmax, color="#777777", lw=0.7, alpha=0.35)
    else:
        ax.set_xlim(0, xmax)
    ax.axvline(1, color="#333333", lw=0.8, ls="--", alpha=0.45)

    for method_i, method in enumerate(METHOD_ORDER):
        method_rows = data[data["method"] == method]
        for _, row in method_rows.iterrows():
            true_x = float(row[metric])
            plot_x = float(row[plot_metric])
            is_clipped = bounded and (true_x < xmin or true_x > xmax)
            y = method_i + tf_offsets[row["TF"]]
            ax.scatter(
                plot_x,
                y,
                s=18 if is_clipped else 15,
                color=b.METHOD_COLORS[method],
                edgecolor="black" if is_clipped else "white",
                linewidth=0.45,
                clip_on=not bounded,
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
                ax.text(text_x, y, text, ha=ha, va="center", fontsize=5.1, color="#111111", clip_on=False, zorder=5)
            else:
                ax.text(
                    plot_x + x_pad,
                    y,
                    row["TF"],
                    ha="left",
                    va="center",
                    fontsize=5.2,
                    color="#222222",
                    clip_on=not bounded,
                    zorder=5,
                )

        vals = method_rows[metric].replace([np.inf, -np.inf], np.nan).dropna()
        if vals.empty:
            continue
        mean_val = float(vals.mean())
        median_val = float(vals.median())
        mean_plot = float(np.clip(mean_val, xmin, xmax)) if bounded else mean_val
        median_plot = float(np.clip(median_val, xmin, xmax)) if bounded else median_val
        summary_rows.extend(
            [
                {
                    "benchmark": benchmark,
                    "method": method,
                    "metric": label,
                    "statistic": "mean",
                    "value": mean_val,
                    "plotted_value": mean_plot,
                    "n_TFs": int(vals.shape[0]),
                },
                {
                    "benchmark": benchmark,
                    "method": method,
                    "metric": label,
                    "statistic": "median",
                    "value": median_val,
                    "plotted_value": median_plot,
                    "n_TFs": int(vals.shape[0]),
                },
            ]
        )
        ax.scatter(
            median_plot,
            method_i,
            marker="D",
            s=42,
            facecolor="white",
            edgecolor="#111111",
            linewidth=0.9,
            clip_on=not bounded,
            zorder=7,
        )
        ax.scatter(
            mean_plot,
            method_i,
            marker="|",
            s=210,
            color="#111111",
            linewidth=1.4,
            clip_on=not bounded,
            zorder=8,
        )

    ax.set_title(label, fontsize=12, pad=8)
    ax.set_xlabel(label)
    ax.set_ylabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_benchmark(data: pd.DataFrame, benchmark: str, bounded: bool, a: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    bench = data[data["benchmark"] == benchmark].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.6 * max(1.0, len(METHOD_ORDER) / 8)), sharey=True)  # Fig5_GRN_v3
    summary_rows: list[dict] = []
    clipped_rows: list[dict] = []
    for ax, (metric, label) in zip(axes, METRICS):
        plot_points_and_stats(ax, bench, benchmark, metric, label, bounded, clipped_rows, summary_rows)

    method_labels = [VIOLIN_METHOD_LABELS[m] for m in METHOD_ORDER]
    axes[0].set_yticks(range(len(method_labels)))
    axes[0].set_yticklabels(method_labels, fontsize=9)
    handles = [
        Line2D([0], [0], marker="D", linestyle="", markerfacecolor="white", markeredgecolor="#111111", markersize=5.5, label="Median"),
        Line2D([0], [0], marker="|", linestyle="", color="#111111", markersize=12, markeredgewidth=1.4, label="Mean"),
    ]
    axes[1].legend(handles=handles, frameon=False, loc="lower right", fontsize=8)
    title_prefix = "bounded " if bounded else ""
    fig.suptitle(f"{BENCHMARK_LABELS[benchmark]} {title_prefix}per-TF performance (sEPR a={a:g})", fontsize=13, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    suffix = "bounded_sEPR_AUPRCr_mean_median_all_methods" if bounded else "sEPR_AUPRCr_mean_median_all_methods"
    save(fig, f"publication_violin_{benchmark}_per_tf_{suffix}{a_suffix(a)}")
    return pd.DataFrame(summary_rows), pd.DataFrame(clipped_rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot sEPR/AUPRCr violin figures.")
    parser.add_argument("--a", type=float, default=A_PRIOR, help="Null-shrinkage pseudocount for sEPR.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = add_sepr(pd.read_csv(METRICS_TABLE), a=args.a)
    sepr_cols = [
        "benchmark",
        "method",
        "method_label",
        "TF",
        "K",
        "M",
        "k",
        "X",
        "p0",
        "lambda",
        "EPR_raw",
        "sEPR",
        "sEPR_a",
        "AUPRC_ratio",
    ]
    sepr_table = suffixed_path(SEPR_TABLE, args.a)
    unbounded_summary_path = suffixed_path(UNBOUNDED_SUMMARY, args.a)
    bounded_summary_path = suffixed_path(BOUNDED_SUMMARY, args.a)
    bounded_clipped_path = suffixed_path(BOUNDED_CLIPPED, args.a)
    data[sepr_cols].to_csv(sepr_table, index=False)

    unbounded_summary = []
    bounded_summary = []
    bounded_clipped = []
    for benchmark in ["TRRUST", "K562_ChIP"]:
        s, _ = plot_benchmark(data, benchmark, bounded=False, a=args.a)
        unbounded_summary.append(s)
        s, c = plot_benchmark(data, benchmark, bounded=True, a=args.a)
        bounded_summary.append(s)
        bounded_clipped.append(c)
    pd.concat(unbounded_summary, ignore_index=True).to_csv(unbounded_summary_path, index=False)
    pd.concat(bounded_summary, ignore_index=True).to_csv(bounded_summary_path, index=False)
    pd.concat(bounded_clipped, ignore_index=True).to_csv(bounded_clipped_path, index=False)
    print(sepr_table)
    print(unbounded_summary_path)
    print(bounded_summary_path)
    print(bounded_clipped_path)


if __name__ == "__main__":
    main()
