#!/usr/bin/env python
"""5b K562 ChIP strip-forests with EPR at a fixed k = 500 instead of k = the TF's number of ChIP targets.

EPR in 5b is precision at k = k_pos (1,435-8,509 ChIP targets per TF), which is larger than the number of targets most
BEELINE runs rank per TF (GRNBOOST2 132-335, the others ~2,300), so for those it measures only which genes are ranked,
not their order. This step re-scores every 5b K562 ChIP (method, TF) at k = 500:

    EPR@500 = (TP in the TF's top 500 targets / 500) / rho,   rho = k_pos / N

Same inputs as the 5b per-TF table (02_reproduce_figures.py): each method's ranked edges
(pub_v2.load_ranked_edges_for_v2, incl. the Yao IRF3 / Replogle FOS supplements), the K562 ChIP reference and the 19,338
valid targets; per TF, edges sorted by EdgeWeight (stable) with duplicates dropped, as in benchmark_grn_systems
.evaluate_one. When a method ranks fewer than 500 targets for a TF, the remaining positions are filled at the expected
hit rate of its unranked targets, (k_pos - TP_emitted) / (N - n_pred), evaluate_one's random-tail rule. The 5b TF sets
are used (requested_tfs: FOS and ELK1 dropped for K562 methods, FOS for THP1). AUPRCr is the 5b value, unchanged.

Writes tables/publication_figures/system_broad_filtered_publication_v2_per_tf_epr_at_500_K562_ChIP.csv and, in
figures/5b/, publication_strip_forest_K562_ChIP_per_tf_EPR500_AUPRCr_all_methods_all_labeled_v2.pdf (+ _summary.csv,
_clipped_points.csv) and ..._v3.pdf (label-placed), drawn by the unchanged 5b plotters (PDF only, as kept for 5b).

Fig5_GRN_v3 addition (EPR fix check): every row of the 5b per-TF table (both benchmarks, every method and TF) is
re-scored here at k = k_pos with the same plain top-k count (epr_at_k: hits[:k]) and must equal the table's TP_at_k /
EPR_raw from benchmark_grn_systems.evaluate_one; the number of rows the pre-fix rule (hits through rank k + 1) changes
is reported. Writes tables/publication_figures/system_broad_filtered_publication_v2_per_tf_epr_check.csv.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
for path in (SCRIPTS_DIR, SCRIPTS_DIR / "capsule_lib" / "01_generate_predictions_and_metrics", SCRIPTS_DIR / "capsule_lib" / "02_visualize_publication_figures"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import os  # noqa: E402

import fig5_paths  # noqa: E402

os.environ.setdefault("MPLCONFIGDIR", str(fig5_paths.OUTPUT_DIR / ".mplconfig"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import benchmark_grn_systems as b  # noqa: E402
import generate_system_broad_filtered_publication_figures_v2 as pub_v2  # noqa: E402

K = 500
BENCHMARK = "K562_ChIP"
METRIC = "EPR_at_500"
TABLE = fig5_paths.PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_epr_at_500_K562_ChIP.csv"
PER_TF = fig5_paths.PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv"
FIG_DIR = fig5_paths.FIGURE_DIR / "5b"
STEM = "publication_strip_forest_K562_ChIP_per_tf_EPR500_AUPRCr_all_methods_all_labeled"
# EPR@500 axis (data range 2026-10-20: 1.22-2.51, no point clipped); AUPRCr keeps 5b's K562 ChIP axis.
EPR500_AXIS = (0.9, 2.7, 1.55)  # xmin, xmax, label threshold (v2 plot_metric)
EPR500_TICKS = np.round(np.arange(1.0, 2.6 + 1e-9, 0.2), 1)


def epr_at_k(edges: pd.DataFrame, positives: set[str], n_universe: int, k: int) -> dict:
    pred = (edges.sort_values("EdgeWeight", ascending=False, kind="mergesort")
            .drop_duplicates(["TF", "target"], keep="first"))
    hits = pred["target"].isin(positives).to_numpy()
    n_pred, k_pos = len(hits), len(positives)
    tp_emitted = int(hits.sum())
    if n_pred >= k:
        tp = float(hits[:k].sum())
    else:
        tail_rate = (k_pos - tp_emitted) / (n_universe - n_pred) if n_universe > n_pred else 0.0
        tp = tp_emitted + (k - n_pred) * tail_rate
    rho = k_pos / n_universe
    return {"k": k, "k_pos": k_pos, "N": n_universe, "rho": rho, "n_pred": n_pred, "TP_emitted": tp_emitted,
            f"TP_at_{k}": tp, METRIC: (tp / k) / rho, f"EPR_max_at_{k}": min(k, k_pos) / k / rho}


def build_table() -> pd.DataFrame:
    import plot_publication_v4_epr_strip_forest as v2

    pub_v2.ALL_METHOD_ORDER = b.METHOD_ORDER  # as 02_reproduce_figures.generate_publication_v2_all_method_tables
    valid = b.load_valid_targets()
    gt = pub_v2.load_reference_ground_truths(valid)[BENCHMARK]
    requested = {tf for m in b.METHOD_ORDER for tf in v2.requested_tfs(BENCHMARK, m)}
    edges = pub_v2.load_ranked_edges_for_v2(valid, requested)
    rows = []
    for method in b.METHOD_ORDER:
        e = edges[method]
        for tf in v2.requested_tfs(BENCHMARK, method):
            positives = set(gt.loc[gt["TF"] == tf, "target"]) & valid
            sub = e[(e["TF"] == tf) & e["target"].isin(valid) & (e["target"] != tf)]
            if not positives or sub.empty:
                continue
            rows.append({"benchmark": BENCHMARK, "method": method, "TF": tf, **epr_at_k(sub, positives, len(valid), K)})
    out = pd.DataFrame(rows)
    pub = pd.read_csv(PER_TF)
    pub = pub[pub["benchmark"].eq(BENCHMARK)][["method", "TF", "k_pos", "n_pred", "EPR", "AUPRC_ratio"]]
    out = out.merge(pub.rename(columns={"k_pos": "k_pos_5b", "n_pred": "n_pred_5b", "EPR": "EPR_at_k_pos"}),
                    on=["method", "TF"], how="left")
    assert out["AUPRC_ratio"].notna().all(), "a (method, TF) is missing from the 5b per-TF table"
    assert (out["k_pos"] == out["k_pos_5b"]).all() and (out["n_pred"] == out["n_pred_5b"]).all(), \
        "positives or ranked-target counts differ from the 5b per-TF table"
    out = out.drop(columns=["k_pos_5b", "n_pred_5b"])
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(TABLE, index=False)
    return out


def plot(table: pd.DataFrame) -> None:
    import seaborn as sns

    import plot_publication_v3_trrust_strip_forest as strip
    import plot_publication_v4_epr_strip_forest as v2
    import plot_strip_forest_v3 as v3

    sns.set_theme(style="white", font_scale=0.9)  # the theme a full 5b run leaves active
    tmp = FIG_DIR / ".epr500_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    xmin, xmax, thr = EPR500_AXIS
    ticks = EPR500_TICKS
    if fig5_paths.SHARED_UNIVERSE:  # 5b_v2: axes for the shared gene universe (also the AUPRCr panel)
        import shared_universe_axes

        shared_universe_axes.apply()
        (xmin, xmax, thr), ticks = shared_universe_axes.EPR500_AXIS, shared_universe_axes.EPR500_TICKS
    auprc = [m for m in v2.BENCHMARK_CONFIGS[BENCHMARK]["metrics"] if m[0] == "AUPRC_ratio"][0]
    v2.BENCHMARK_CONFIGS[BENCHMARK]["metrics"] = [(METRIC, f"EPR (k = {K})", xmin, xmax, thr), auprc]
    v2.X_TICKS[(BENCHMARK, METRIC)] = ticks
    v2.OUT_DIR = v3.OUT_DIR = strip.OUT_DIR = tmp
    strip.ordered_tfs_for_benchmark = v2.ordered_tfs_for_benchmark
    data = table[["benchmark", "method", "TF", METRIC, "AUPRC_ratio"]]
    v2.plot_benchmark(data, BENCHMARK)
    v3.plot_benchmark(data, BENCHMARK)
    old = f"publication_strip_forest_{BENCHMARK}_per_tf_EPR_AUPRCr_all_methods_all_labeled"
    for f in sorted(tmp.iterdir()):
        if f.suffix == ".png":  # 5b keeps the v2 strip-forests as PDF only
            f.unlink()
            continue
        dest = FIG_DIR / f.name.replace(old, STEM)
        shutil.move(str(f), dest)
        print("wrote", dest)
    tmp.rmdir()


def check_epr_at_k_pos() -> pd.DataFrame:
    """Independent recount of the 5b per-TF table's TP_at_k (k = k_pos) for every (benchmark, method, TF)."""
    pub_v2.ALL_METHOD_ORDER = b.METHOD_ORDER
    valid = b.load_valid_targets()
    gts = pub_v2.load_reference_ground_truths(valid)
    per_tf = pd.read_csv(PER_TF)
    edges = pub_v2.load_ranked_edges_for_v2(valid, set(per_tf["TF"]))
    rows = []
    for r in per_tf.itertuples(index=False):
        positives = set(gts[r.benchmark].loc[lambda g: g["TF"] == r.TF, "target"]) & valid
        e = edges[r.method]
        sub = e[(e["TF"] == r.TF) & e["target"].isin(valid) & (e["target"] != r.TF)]
        k = int(r.k_pos)
        new = epr_at_k(sub, positives, len(valid), k)
        hits = (sub.sort_values("EdgeWeight", ascending=False, kind="mergesort")
                .drop_duplicates(["TF", "target"], keep="first")["target"].isin(positives).to_numpy())
        old_tp = float(hits[:k + 1].sum()) if len(hits) > k else np.nan  # pre-fix rule (only differs when n_pred > k)
        rows.append({"benchmark": r.benchmark, "method": r.method, "TF": r.TF, "k_pos": k, "n_pred": len(hits),
                     "TP_at_k_table": r.TP_at_k, "TP_at_k_recount": new[f"TP_at_{k}"], "EPR_table": r.EPR_raw,
                     "EPR_recount": new[METRIC], "TP_at_k_prefix_rule": old_tp})
    out = pd.DataFrame(rows)
    out.to_csv(fig5_paths.PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_epr_check.csv", index=False)
    d = (out["TP_at_k_table"] - out["TP_at_k_recount"]).abs().max()
    assert len(out) == len(per_tf) and d < 1e-9 and np.allclose(out["EPR_table"], out["EPR_recount"]), \
        f"evaluate_one TP_at_k differs from the plain top-k recount (max |delta| {d})"
    changed = (out["TP_at_k_prefix_rule"] - out["TP_at_k_recount"]).abs() > 1e-9
    print(f"EPR check: {len(out)} per-TF rows, evaluate_one TP_at_k == plain top-k recount (max |delta| {d:.3g}); "
          f"the pre-fix rule (hits through rank k + 1) would change {int(changed.sum())} of them")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--table-only", action="store_true", help="compute and print the EPR@500 table; no figures")
    args = parser.parse_args()
    check_epr_at_k_pos()
    table = build_table()
    with pd.option_context("display.width", 200, "display.max_rows", 100):
        print(table[["method", "TF", "k_pos", "n_pred", "TP_emitted", f"TP_at_{K}", METRIC, f"EPR_max_at_{K}",
                     "EPR_at_k_pos", "AUPRC_ratio"]].round(3).to_string(index=False))
        print(f"\n{METRIC} range {table[METRIC].min():.3f}-{table[METRIC].max():.3f}; median by method:")
        print(table.groupby("method", sort=False)[METRIC].median().round(2).to_string())
    if not args.table_only:
        plot(table)
    print("EPR500_OK")


if __name__ == "__main__":
    main()
