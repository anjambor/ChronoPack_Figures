# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import benchmark_grn_systems as b
from generate_system_broad_filtered_publication_figures_v2 import (
    ALL_METHOD_ORDER,
    BASE_TFS_BY_BENCHMARK,
    SCENARIO,
    THP1_EXTRA_TFS,
    load_ranked_edges_for_v2,
    load_reference_ground_truths,
)
from plot_publication_v2_violin_sepr import add_sepr


BASE_DIR = fig5_paths.BENCH_ROOT
FIG_BASE = BASE_DIR / "benchmark_outputs" / "figures"
V3_DIR = fig5_paths.CAPSULE_INPUT_DIR / "direct_figure_inputs" / "publication_figures_v3"
V4_DIR = FIG_BASE
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables"
PUBLICATION_TABLE_DIR = TABLE_DIR / "publication_figures"
V2_TABLE_DIR = PUBLICATION_TABLE_DIR
V4_TABLE_DIR = PUBLICATION_TABLE_DIR

K562_DROP_TFS = {"FOS", "ELK1"}
THP1_DROP_TFS = {"FOS"}
K562_METHODS_ALL = {
    "ChronoPack_K562",
    "GENIE3",
    "SINCERITIES",
    "GRNBOOST2",
    "PPCOR",
    "Replogle2022_K562",
    *fig5_paths.PRIOR_METHODS,
}
THP1_METHODS_ALL = {"ChronoPack_THP1", "Yao2023_THP1"}


def method_context(method: str) -> str:
    return "THP1" if method in THP1_METHODS_ALL else "K562"


def v4_requested_tfs(benchmark: str, method: str) -> list[str]:
    tfs = list(BASE_TFS_BY_BENCHMARK[benchmark])
    if method_context(method) == "THP1":
        tfs.extend(THP1_EXTRA_TFS)
        drop = THP1_DROP_TFS
    else:
        drop = K562_DROP_TFS
    return [tf for tf in dict.fromkeys(tfs) if tf not in drop]


def v4_ordered_tfs_for_benchmark(benchmark: str) -> list[str]:
    out: list[str] = []
    for method in b.METHOD_ORDER:
        out.extend(v4_requested_tfs(benchmark, method))
    return list(dict.fromkeys(out))


def filter_metric_rows(data: pd.DataFrame) -> pd.DataFrame:
    keep = [
        tf in v4_requested_tfs(benchmark, method)
        for benchmark, method, tf in zip(data["benchmark"], data["method"], data["TF"])
    ]
    return data.loc[keep].copy()


def filter_curve_rows(curves: pd.DataFrame) -> pd.DataFrame:
    keep = [
        tf in v4_requested_tfs(benchmark, method)
        for benchmark, method, tf in zip(curves["benchmark"], curves["method"], curves["TF"])
    ]
    return curves.loc[keep].copy()


def generate_6b() -> None:
    import plot_publication_v3_trrust_strip_forest as strip
    import plot_publication_v3_trrust_violin_sepr as violin

    out_dir = V4_DIR / "5b"
    metrics = add_sepr(pd.read_csv(strip.METRICS_TABLE), a=strip.A_PRIOR)
    metrics = filter_metric_rows(metrics)

    strip.OUT_DIR = out_dir
    strip.ordered_tfs_for_benchmark = v4_ordered_tfs_for_benchmark
    for benchmark in ["TRRUST", "K562_ChIP"]:
        strip.plot_benchmark(metrics, benchmark)

    violin.OUT_DIR = out_dir
    violin.ordered_tfs_for_benchmark = v4_ordered_tfs_for_benchmark
    for benchmark in ["TRRUST", "K562_ChIP"]:
        violin.plot_benchmark(metrics, benchmark)


def generate_6c_macro() -> None:
    import plot_publication_v3_mean_macro_pr as macro

    macro.OUT_DIR = V4_DIR / "5c"
    curves = pd.read_csv(macro.CURVES_TABLE)
    curves = filter_curve_rows(curves)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        macro.plot_macro_pr_for_benchmark(curves, benchmark)


def build_v4_micro_curves() -> tuple[pd.DataFrame, pd.DataFrame]:
    valid_targets = b.load_valid_targets()
    gts = load_reference_ground_truths(valid_targets)
    requested_tfs = set()
    for benchmark in ["TRRUST", "K562_ChIP"]:
        for method in ALL_METHOD_ORDER:
            requested_tfs.update(v4_requested_tfs(benchmark, method))

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
                for tf in v4_requested_tfs(benchmark, method)
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
                    "method_label": b.METHOD_LABELS[method],
                    "method_context": method_context(method),
                    "n_tfs_eval": len(tf_set),
                    "tfs": ",".join([tf for tf in v4_requested_tfs(benchmark, method) if tf in tf_set]),
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
    V4_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(V4_TABLE_DIR / "system_broad_filtered_publication_v4_micro_pr_matched_metrics.csv", index=False)
    curves.to_csv(V4_TABLE_DIR / "system_broad_filtered_publication_v4_micro_pr_matched_curves.csv", index=False)
    return metrics, curves


def generate_6c_micro() -> None:
    import plot_publication_v3_micro_pr as micro

    metrics_path = V4_TABLE_DIR / "system_broad_filtered_publication_v4_micro_pr_matched_metrics.csv"
    curves_path = V4_TABLE_DIR / "system_broad_filtered_publication_v4_micro_pr_matched_curves.csv"
    if metrics_path.exists() and curves_path.exists():
        metrics = pd.read_csv(metrics_path)
        curves = pd.read_csv(curves_path)
    else:
        metrics, curves = build_v4_micro_curves()

    micro.OUT_DIR = V4_DIR / "5c"
    print(metrics[["benchmark", "method", "n_tfs_eval", "tfs", "k_pos", "TP_emitted", "AUPRC_ratio"]].to_string(index=False))
    for benchmark in ["TRRUST", "K562_ChIP"]:
        micro.plot_micro_pr_for_benchmark(curves, benchmark)


def generate_6d() -> None:
    import plot_publication_v3_6d_trrust_pr_curves as pr6d
    import plot_publication_v3_trrust_core_heatmap_reference_palette as heat6d

    pr6d.OUT_DIR = V4_DIR / "5d"
    curves = pd.read_csv(pr6d.TABLE)
    for tf in ["RELA", "NFKB1", "JUN"]:
        pr6d.plot_tf(curves, tf, ymax=pr6d.YLIMS[tf])
        pr6d.plot_tf(curves, tf, ymax=pr6d.VARIANT_YLIM, suffix="_ymax1p0")

    heat6d.OUT_DIR = V4_DIR / "5d"
    heat6d.SOURCE_OUT = V4_TABLE_DIR / "system_broad_filtered_publication_v4_trrust_core_tf_metrics_reference_palette.csv"
    heat6d.main()

    for path in (V3_DIR / "6d").glob("expr-curve_*.pdf"):
        (V4_DIR / "5d").mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, V4_DIR / "5d" / path.name)


def generate_6e() -> None:
    import plot_regulator_prioritization_v2 as reg
    import plot_regulator_prioritization_v2_composite_stacked as comp
    import plot_regulator_prioritization_v2_nes_dot_comparison as nes

    k562_eval = [tf for tf in reg.K562_EVAL_TFS if tf not in K562_DROP_TFS]
    thp1_eval = [tf for tf in reg.THP1_EVAL_TFS if tf not in THP1_DROP_TFS]

    def requested_tfs(method: str) -> list[str]:
        return thp1_eval if method.endswith("_THP1") else k562_eval

    def all_requested_tfs() -> list[str]:
        return list(dict.fromkeys(k562_eval + thp1_eval))

    def load_composition_pivot() -> pd.DataFrame:
        data = pd.read_csv(comp.TOPK_PATH)
        sub = data[data["top_k"].eq(10000)].copy()
        keep = [tf in requested_tfs(method) for method, tf in zip(sub["method"], sub["TF"])]
        sub = sub.loc[keep].copy()
        plot = sub.groupby(["method", "TF"], as_index=False).agg(edge_share=("edge_share", "sum"))
        comparator_totals = (
            plot[plot["method"].isin(comp.PANEL_B_COMPARATOR_METHODS)]
            .groupby("method")["edge_share"]
            .sum()
            .sort_values(ascending=False)
        )
        method_order = comp.PANEL_B_CHRONO_METHODS + comparator_totals.index.tolist()
        return (
            plot.pivot_table(index="method", columns="TF", values="edge_share", aggfunc="sum", fill_value=0)
            .reindex(index=method_order, columns=all_requested_tfs(), fill_value=0)
        )

    nes.K562_TFS = k562_eval
    comp.FIG_DIR = V4_DIR / "5e"
    comp.SUMMARY_OUT = V4_TABLE_DIR / "regulator_prioritization_v4_chrono_vs_best_comparator_nes.csv"
    comp.requested_tfs = requested_tfs
    comp.all_requested_tfs = all_requested_tfs
    comp.load_composition_pivot = load_composition_pivot
    comp.main()


def main() -> None:
    V4_DIR.mkdir(parents=True, exist_ok=True)
    V4_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    # Fig5_GRN_v3: generate_6b (sEPR strip-forests and sEPR violins) dropped at the user's request (2026-10-22)
    generate_6c_macro()
    build_v4_micro_curves()  # Fig5_GRN_v3: micro-PR table only; the 5c micro-PR figures dropped (user request, 2026-10-22)
    generate_6d()
    generate_6e()
    print("Wrote publication v4 figures to", V4_DIR)


if __name__ == "__main__":
    main()
