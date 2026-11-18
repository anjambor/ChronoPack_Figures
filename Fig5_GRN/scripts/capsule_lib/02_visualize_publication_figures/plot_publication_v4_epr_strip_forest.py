# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd

import benchmark_grn_systems as b
import plot_publication_v3_trrust_strip_forest as strip


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables" / "publication_figures"
OUT_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "5b"
METRICS_TABLE = TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv"

K562_DROP_TFS = {"FOS", "ELK1"}
THP1_DROP_TFS = {"FOS"}
THP1_METHODS = {"ChronoPack_THP1", "Yao2023_THP1"}
BASE_TFS_BY_BENCHMARK = {
    "TRRUST": ["RELA", "NFKB1", "JUN", "ELK1", "CREB1", "ATF1", "FOS", "ATF2"],
    "K562_ChIP": ["RELA", "JUN", "ELK1", "CREB1", "ATF1", "FOS", "ATF2"],
}
THP1_EXTRA_TFS = ["REL", "IRF3"]
# Fig5 addition: x tick labels per (benchmark, metric); None keeps the strip-plot default (integers).
# Tick marks are drawn on every panel (the active seaborn "white" theme hides them otherwise).
X_TICKS = {
    ("K562_ChIP", "EPR"): np.round(np.arange(1.0, 2.2 + 1e-9, 0.2), 1),
    ("K562_ChIP", "AUPRC_ratio"): np.round(np.arange(1.0, 2.0 + 1e-9, 0.2), 1),
}
BENCHMARK_CONFIGS = {
    "TRRUST": {
        "label": "TRRUST",
        "metrics": [
            ("EPR", "EPR", -0.5, 10.0, 2.0),
            ("AUPRC_ratio", "AUPRCr", 0.5, 5.5, 1.4),
        ],
    },
    "K562_ChIP": {
        "label": "K562 ChIP",
        "metrics": [
            ("EPR", "EPR", 0.9, 2.3, 1.55),
            ("AUPRC_ratio", "AUPRCr", 0.9, 2.0, 1.45),
        ],
    },
}


def method_context(method: str) -> str:
    return "THP1" if method in THP1_METHODS else "K562"


def requested_tfs(benchmark: str, method: str) -> list[str]:
    tfs = list(BASE_TFS_BY_BENCHMARK[benchmark])
    if method_context(method) == "THP1":
        tfs.extend(THP1_EXTRA_TFS)
        drop = THP1_DROP_TFS
    else:
        drop = K562_DROP_TFS
    return [tf for tf in dict.fromkeys(tfs) if tf not in drop]


def ordered_tfs_for_benchmark(benchmark: str) -> list[str]:
    out: list[str] = []
    for method in b.METHOD_ORDER:
        out.extend(requested_tfs(benchmark, method))
    return list(dict.fromkeys(out))


def filter_rows(data: pd.DataFrame) -> pd.DataFrame:
    keep = [
        tf in requested_tfs(benchmark, method)
        for benchmark, method, tf in zip(data["benchmark"], data["method"], data["TF"])
    ]
    return data.loc[keep].copy()


def save(fig: plt.Figure, out_stem: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{out_stem}.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{out_stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_benchmark(all_data: pd.DataFrame, benchmark: str) -> None:
    config = BENCHMARK_CONFIGS[benchmark]
    data = all_data[all_data["benchmark"].eq(benchmark) & all_data["method"].isin(b.METHOD_ORDER)].copy()
    y_values, fig_height, (k562_label_y, thp1_label_y) = b.strip_rows(b.METHOD_ORDER, 4.65)  # Fig5_GRN_v3
    y_positions = dict(zip(b.METHOD_ORDER, y_values))
    y_ticks = [y_positions[m] for m in b.METHOD_ORDER]
    y_labels = [strip.method_label(m) for m in b.METHOD_ORDER]

    fig, axes = plt.subplots(1, 2, figsize=(10.6, fig_height), sharey=True)
    clipped_rows: list[dict] = []
    summary_rows: list[dict] = []
    for ax, (metric, title, xmin, xmax, label_threshold) in zip(axes, config["metrics"]):
        strip.plot_metric(
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
        ticks = X_TICKS.get((benchmark, metric))
        if ticks is not None:
            ax.set_xticks(ticks)
            if not all(float(t).is_integer() for t in ticks):  # Fig5: integer ticks (shared-universe TRRUST) keep plain labels
                ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.tick_params(axis="x", bottom=True, length=3.5, width=0.8, color=ax.spines["bottom"].get_edgecolor())

    for ax in axes:
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_labels, fontsize=11.0)
        ax.set_ylim(max(y_values) + 0.55, -0.62)
    axes[0].text(-0.30, k562_label_y, "K562", transform=axes[0].transAxes, ha="center", va="center", rotation=90, fontsize=9.8, color="#555555")
    axes[0].text(-0.30, thp1_label_y, "THP1", transform=axes[0].transAxes, ha="center", va="center", rotation=90, fontsize=9.8, color="#555555")

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

    out_stem = f"publication_strip_forest_{benchmark}_per_tf_EPR_AUPRCr_all_methods_all_labeled_v2"
    save(fig, out_stem)
    pd.DataFrame(clipped_rows).to_csv(OUT_DIR / f"{out_stem}_clipped_points.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(OUT_DIR / f"{out_stem}_summary.csv", index=False)
    print(OUT_DIR / f"{out_stem}.pdf")
    print(OUT_DIR / f"{out_stem}.png")
    if clipped_rows:
        print(pd.DataFrame(clipped_rows).to_string(index=False))


def main() -> None:
    strip.OUT_DIR = OUT_DIR
    strip.ordered_tfs_for_benchmark = ordered_tfs_for_benchmark
    data = filter_rows(pd.read_csv(METRICS_TABLE))
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_benchmark(data, benchmark)


if __name__ == "__main__":
    main()
