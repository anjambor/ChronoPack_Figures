# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import benchmark_grn_systems as b
from generate_system_broad_filtered_publication_figures_v2 import (
    ALL_METHOD_ORDER,
    BASE_TFS_BY_BENCHMARK,
    METHOD_CONTEXT,
    SCENARIO,
    THP1_EXTRA_TFS,
    load_ranked_edges_for_v2,
    load_reference_ground_truths,
    method_requested_tfs,
)


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables"
OUT_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "5c"
V3_TABLE_DIR = TABLE_DIR / "publication_figures"
MATCHED_METRICS_TABLE = V3_TABLE_DIR / "system_broad_filtered_publication_v3_micro_pr_matched_metrics.csv"
MATCHED_CURVES_TABLE = V3_TABLE_DIR / "system_broad_filtered_publication_v3_micro_pr_matched_curves.csv"

BENCHMARK_LABELS = {
    "TRRUST": "TRRUST",
    "K562_ChIP": "K562 ChIP",
}
METHOD_COLORS = {
    "ChronoPack_K562": "#D99000",
    "GENIE3": "#3D8B37",
    "Replogle2022_K562": "#CC79A7",
    "ChronoPack_THP1": "#D55E00",
    "Yao2023_THP1": "#3366AA",
    **fig5_paths.PRIOR_COLORS,
}
METHOD_LABELS = {
    "ChronoPack_K562": "ChronoPack | K562 TNF",
    "GENIE3": "GENIE3 | K562 TNF",
    "Replogle2022_K562": "Replogle2022 | K562 perturb",
    "ChronoPack_THP1": "ChronoPack | THP1 LPS",
    "Yao2023_THP1": "Yao2023 | THP1 perturb",
    **{m: f"{label} | K562 TNF" for m, label in fig5_paths.PRIOR_LABELS.items()},
}
YLIMS = {
    "TRRUST": (0.0, 0.35),
    "K562_ChIP": (0.0, 0.5),
}
TITLE_SIZE = 16
LABEL_SIZE = 15
TICK_SIZE = 15


def save(fig: plt.Figure, stem: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def display_curve(d: pd.DataFrame, max_points: int = 1200) -> pd.DataFrame:
    d = d.sort_values("k").dropna(subset=["recall", "precision"])
    if len(d) <= max_points:
        return d
    idx = np.linspace(0, len(d) - 1, max_points).round().astype(int)
    idx = np.unique(np.r_[0, idx, len(d) - 1])
    return d.iloc[idx]


def build_matched_micro_curves() -> tuple[pd.DataFrame, pd.DataFrame]:
    valid_targets = b.load_valid_targets()
    gts = load_reference_ground_truths(valid_targets)
    requested_tfs = set()
    for benchmark in ["TRRUST", "K562_ChIP"]:
        for method in ALL_METHOD_ORDER:
            requested_tfs.update(method_requested_tfs(benchmark, method))

    method_edges = load_ranked_edges_for_v2(valid_targets, requested_tfs)
    metric_rows = []
    curve_rows = []
    for benchmark, gt in gts.items():
        gt_tfs = set(gt["TF"])
        for method in ALL_METHOD_ORDER:
            edges = method_edges[method]
            available_tfs = set(edges["TF"])
            tf_set = {
                tf
                for tf in method_requested_tfs(benchmark, method)
                if tf in gt_tfs and tf in available_tfs
            }
            if not tf_set:
                continue
            row, curve = b.evaluate_one(edges, gt, tf_set, valid_targets, curve_points=500)
            row.update(
                {
                    "scenario": SCENARIO,
                    "benchmark": benchmark,
                    "scope": "network",
                    "TF": "ALL",
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "method_context": METHOD_CONTEXT[method],
                    "n_tfs_eval": len(tf_set),
                    "tfs": ",".join([tf for tf in method_requested_tfs(benchmark, method) if tf in tf_set]),
                }
            )
            metric_rows.append(row)
            c = curve.copy()
            c["scenario"] = SCENARIO
            c["benchmark"] = benchmark
            c["scope"] = "network"
            c["TF"] = "ALL"
            c["method"] = method
            curve_rows.append(c)

    metrics = pd.DataFrame(metric_rows)
    curves = pd.concat(curve_rows, ignore_index=True)
    V3_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(MATCHED_METRICS_TABLE, index=False)
    curves.to_csv(MATCHED_CURVES_TABLE, index=False)
    return metrics, curves


def load_or_build_matched_micro_curves() -> tuple[pd.DataFrame, pd.DataFrame]:
    if MATCHED_METRICS_TABLE.exists() and MATCHED_CURVES_TABLE.exists():
        return pd.read_csv(MATCHED_METRICS_TABLE), pd.read_csv(MATCHED_CURVES_TABLE)
    return build_matched_micro_curves()


def plot_micro_pr_for_benchmark(curves: pd.DataFrame, benchmark: str) -> None:
    panel = curves[
        curves["scenario"].eq(SCENARIO)
        & curves["benchmark"].eq(benchmark)
        & curves["scope"].eq("network")
        & curves["TF"].eq("ALL")
        & curves["method"].isin(ALL_METHOD_ORDER)
    ].copy()

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for method in ALL_METHOD_ORDER:
        d = panel[panel["method"].eq(method)].sort_values("k")
        if d.empty:
            continue
        observed = display_curve(d[~d["is_expectation"]])
        expected = display_curve(d[d["is_expectation"]], max_points=500)
        color = METHOD_COLORS[method]
        label = METHOD_LABELS[method]
        if len(observed):
            ax.plot(
                observed["recall"],
                observed["precision"],
                drawstyle="steps-post",
                lw=2.4,
                color=color,
                label=label,
            )
        if len(expected):
            tail = pd.concat([observed.tail(1), expected]).drop_duplicates(["recall", "precision"])
            ax.plot(tail["recall"], tail["precision"], lw=2.0, ls="--", color=color, alpha=0.72)

    ax.set_xlim(0, 0.25 if benchmark == "TRRUST" else 1)
    ax.set_ylim(*YLIMS[benchmark])
    if benchmark == "TRRUST":
        ax.set_yticks([0.0, 0.1, 0.2, 0.3])
    ax.set_box_aspect(1)
    ax.set_xlabel("Recall", fontsize=LABEL_SIZE)
    ax.set_ylabel("Precision", fontsize=LABEL_SIZE)
    ax.set_title(BENCHMARK_LABELS[benchmark], fontsize=TITLE_SIZE, pad=10)
    ax.tick_params(axis="both", labelsize=TICK_SIZE)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if benchmark == "TRRUST":
        ax.legend(frameon=False, fontsize=6.5, loc="upper right", handlelength=1.6, borderaxespad=0.2)
    fig.tight_layout()
    save(fig, f"publication_micro_pr_{benchmark}")
    print(OUT_DIR / f"publication_micro_pr_{benchmark}.pdf")
    print(OUT_DIR / f"publication_micro_pr_{benchmark}.png")


def main() -> None:
    metrics, curves = load_or_build_matched_micro_curves()
    print(metrics[["benchmark", "method", "n_tfs_eval", "tfs", "k_pos", "TP_emitted", "AUPRC_ratio"]].to_string(index=False))
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_micro_pr_for_benchmark(curves, benchmark)


if __name__ == "__main__":
    main()
