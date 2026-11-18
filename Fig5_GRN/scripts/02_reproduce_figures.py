#!/usr/bin/env python
"""Rebuild the GRN benchmark tables and publication-v4 figures (panels 5b-5e) for Fig5.

Fig5_GRN_v4: as Fig5_GRN_v3, with step 01's new ChronoPack TF rankings. Run once per method set (FIG5_METHOD_SET=expression, the default, and with_prior); with_prior adds
the motif-prior methods to every cross-method panel and writes to output/with_prior/benchmark_outputs/.

Port of grn_inference_v2/reproduction/scripts/02_visualize_publication_figures/
reproduce_publication_figures_v4.py. Differences from the capsule runner:

- ChronoPack_K562 / ChronoPack_THP1 predictions are read from the chronopack package
  output of 01_run_chronopack_grn.py (ChronoPack r9 data) instead of the capsule's
  re-implementation (explore_thp1_lfdr_recalibration.R + load_chronopack_beta_lfdr).
  The package applies the same edge weighting, normalize_edges and break_ties.
- Comparator predictions (BEELINE, Replogle, Yao) and ground truths are read from the
  copies in input/capsule/external/ via the capsule loaders, unchanged.
- Every stage is regenerated (no skip-if-exists); output/benchmark_outputs/ is cleared first.
- Adds the FOS/ELK1-filtered K562 ChIP macro-PR zooms (x 0-0.5, y 0.2-0.5, previously made ad
  hoc, and y 0.2-0.6).

``--zoom-only`` regenerates only those zooms, from the existing PR-curve table, and ``--strip-only`` only the
5b EPR/AUPRCr strip-forests (v2 and the label-placed v3), from the existing per-TF metrics table
(nothing is cleared).
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
for path in (SCRIPTS_DIR, SCRIPTS_DIR / "capsule_lib" / "01_generate_predictions_and_metrics", SCRIPTS_DIR / "capsule_lib" / "02_visualize_publication_figures"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import fig5_paths
import shared_universe_axes

os.environ.setdefault("MPLCONFIGDIR", str(fig5_paths.OUTPUT_DIR / ".mplconfig"))

import matplotlib
import pandas as pd

TABLE_DIR = fig5_paths.TABLE_DIR
PUBLICATION_TABLE_DIR = fig5_paths.PUBLICATION_TABLE_DIR
V4_FIG_DIR = fig5_paths.FIGURE_DIR
RANKED_EDGES = TABLE_DIR / "regulator_prioritization_method_specific_ranked_edges.csv"


def load_chronopack_package_edges(name: str, context: str, valid_targets: set[str]) -> pd.DataFrame:
    """Package GRN edges in the capsule's (TF, target, EdgeWeight, method) schema, in package rank order."""
    edges = pd.read_csv(fig5_paths.CHRONOPACK_OUTPUT_DIR / context / "grn.edges.csv")
    if set(edges["method"]) != {name}:
        raise ValueError(f"Unexpected method labels in {context} edges: {sorted(set(edges['method']))}")
    edges = edges.sort_values("rank", kind="mergesort")
    edges = edges.rename(columns={"regulator": "TF", "score": "EdgeWeight"})
    edges = edges[edges["target"].astype(str).isin(valid_targets)]
    return edges[["TF", "target", "EdgeWeight", "method"]].reset_index(drop=True)


def install_chronopack_loader(compare) -> None:
    c = compare.c

    def load_methods(valid_targets: set[str]) -> dict[str, pd.DataFrame]:
        methods = {name: c.load_beeline_method(name, valid_targets) for name in ["GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"]}
        for context, (name, _source, _is_ligand) in fig5_paths.CONTEXTS.items():
            methods[name] = load_chronopack_package_edges(name, context, valid_targets)
        methods["Replogle2022_K562"] = c.load_replogle(valid_targets)
        methods["Yao2023_THP1"] = c.load_yao(valid_targets)
        if fig5_paths.WITH_PRIOR:  # Fig5_GRN_v3: motif-prior networks (../Run_PriorGRN)
            methods.update({name: c.load_prior_method(name, valid_targets) for name in fig5_paths.PRIOR_METHODS})
        return {method: methods[method] for method in compare.METHOD_ORDER if method in methods}

    compare.load_recalibrated_methods = load_methods


def generate_regulator_prioritization_base_tables(compare) -> None:
    valid_targets = compare.c.load_valid_targets()
    reference_tfs = compare.c.load_reference_tfs(valid_targets)
    methods = compare.load_recalibrated_methods(valid_targets)
    universes = compare.method_universes(methods, reference_tfs)
    rank_rows = []
    for method in compare.METHOD_ORDER:
        universe_tfs, universe_targets = universes[method]
        pred = methods[method]
        pred = pred[pred["TF"].isin(universe_tfs) & pred["target"].isin(universe_targets)].copy()
        pred = (
            pred.sort_values("EdgeWeight", ascending=False, kind="mergesort")
            .drop_duplicates(["TF", "target"], keep="first")
            .reset_index(drop=True)
        )
        pred["rank"] = range(1, len(pred) + 1)
        pred["percentile"] = 1 - (pred["rank"] - 1) / max(len(pred), 1)
        rank_rows.append(pred[["method", "TF", "target", "EdgeWeight", "rank", "percentile"]])
        print(f"  ranked edges {method}: {len(pred)} edges, {pred['TF'].nunique()} TFs")
    pd.concat(rank_rows, ignore_index=True).to_csv(RANKED_EDGES, index=False)


def generate_regulator_prioritization_v2_tables(reg_v2) -> None:
    reg_v2.load_topk_counts().to_csv(TABLE_DIR / "regulator_prioritization_v2_fixed_topk_tf_proportions_full.csv", index=False)
    metrics, _ = reg_v2.compute_rank_enrichment(include_curves=False)
    metrics.to_csv(TABLE_DIR / "regulator_prioritization_v2_tf_rank_enrichment.csv", index=False)


def generate_publication_v2_base_tables(pub_v2, b) -> None:
    pub_v2.V2_TABLE_DIR = PUBLICATION_TABLE_DIR
    valid_targets = b.load_valid_targets()
    gts = pub_v2.load_reference_ground_truths(valid_targets)
    requested = set()
    for benchmark in ["TRRUST", "K562_ChIP"]:
        for method in pub_v2.ALL_METHOD_ORDER:
            requested.update(pub_v2.method_requested_tfs(benchmark, method))
    method_edges = pub_v2.load_ranked_edges_for_v2(valid_targets, requested)
    _, curves, _ = pub_v2.build_v2_metrics(method_edges, gts, valid_targets)
    curves.to_csv(PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_pr_curves.csv", index=False)


def generate_publication_v2_all_method_tables(pub_v2, b) -> None:
    pub_v2.V2_TABLE_DIR = PUBLICATION_TABLE_DIR
    pub_v2.ALL_METHOD_ORDER = b.METHOD_ORDER
    pub_v2.METHOD_CONTEXT = {
        "ChronoPack_K562": "K562",
        "GENIE3": "K562",
        "SINCERITIES": "K562",
        "GRNBOOST2": "K562",
        "PPCOR": "K562",
        "Replogle2022_K562": "K562",
        "ChronoPack_THP1": "THP1",
        "Yao2023_THP1": "THP1",
        **{m: "K562" for m in fig5_paths.PRIOR_METHODS},  # Fig5_GRN_v3
    }
    valid_targets = b.load_valid_targets()
    gts = pub_v2.load_reference_ground_truths(valid_targets)
    requested = set()
    for benchmark in ["TRRUST", "K562_ChIP"]:
        for method in b.METHOD_ORDER:
            requested.update(pub_v2.method_requested_tfs(benchmark, method))
    method_edges = pub_v2.load_ranked_edges_for_v2(valid_targets, requested)
    metrics, _, _ = pub_v2.build_v2_metrics(method_edges, gts, valid_targets)
    metrics.to_csv(PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv", index=False)


def generate_atf2_heatmap_variant() -> None:
    import plot_publication_v3_trrust_core_heatmap_reference_palette as heat

    heat.TABLE_DIR = PUBLICATION_TABLE_DIR
    heat.METRICS_TABLE = PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv"
    heat.OUT_DIR = V4_FIG_DIR / "5d"
    heat.CORE_TFS = ["RELA", "NFKB1", "JUN", "ATF2"]
    heat.OUT_STEM = "publication_heatmap_TRRUST_RELA_NFKB1_JUN_ATF2_EPR_AUPRCr"
    heat.SOURCE_OUT = PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v4_trrust_core_tf_metrics_reference_palette_RELA_NFKB1_JUN_ATF2.csv"
    heat.FIGSIZE = (9.6, 4.0)
    heat.main()


# K562 ChIP macro-PR zooms (x 0-0.5): output suffix -> (y limits, y ticks, full-panel style).
# y ticks None = the plot's default rule (0.2 / 0.4 within the limits). Full-panel style = matplotlib
# defaults, as active when v4 draws publication_mean_macro_pr_K562_ChIP (tick marks, black text);
# otherwise the zoom inherits the seaborn "white" theme left by the ATF2 heatmap (no tick marks).
K562_MACRO_ZOOMS = {
    "_xmax0p5_ymin0p2": ((0.2, 0.5), None, False),  # formerly an ad hoc scratchpad script
    "_xmax0p5_ymin0p2_ymax0p6": ((0.2, 0.6), [0.2, 0.4, 0.6], True),
}


def generate_k562_macro_zoom(v4) -> None:
    """FOS/ELK1-filtered K562 ChIP macro-PR zooms (K562_MACRO_ZOOMS; shared_universe_axes' for FIG5_UNIVERSE=shared)."""
    import plot_publication_v3_mean_macro_pr as macro

    macro.OUT_DIR = V4_FIG_DIR / "5c"
    curves = v4.filter_curve_rows(pd.read_csv(macro.CURVES_TABLE))
    zooms = shared_universe_axes.K562_MACRO_ZOOMS if fig5_paths.SHARED_UNIVERSE else K562_MACRO_ZOOMS
    for suffix, (ylim, yticks, full_panel_style) in zooms.items():
        with matplotlib.rc_context(matplotlib.rcParamsDefault if full_panel_style else None):
            macro.plot_macro_pr_for_benchmark(curves, "K562_ChIP", xlim=(0.0, 0.5), ylim=ylim, yticks=yticks, out_suffix=suffix)


def generate_strip_forest_v3() -> None:
    """5b strip-forests with clash-avoiding labels (scripts/plot_strip_forest_v3.py); PDF only."""
    import plot_strip_forest_v3 as v3

    v3.OUT_DIR = V4_FIG_DIR / "5b"
    v3.main()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--zoom-only", action="store_true", help="only regenerate the 5c K562 ChIP macro-PR zooms")
    parser.add_argument("--strip-only", action="store_true", help="only regenerate the 5b EPR/AUPRCr strip-forests")
    args = parser.parse_args()
    if args.strip_only:
        import plot_publication_v4_epr_strip_forest as epr_strip
        import seaborn as sns

        if fig5_paths.SHARED_UNIVERSE:
            shared_universe_axes.apply()

        # In a full run the strip-forests follow v4.main(), which leaves this theme active; match it.
        sns.set_theme(style="white", font_scale=0.9)
        epr_strip.OUT_DIR = V4_FIG_DIR / "5b"
        epr_strip.TABLE_DIR = PUBLICATION_TABLE_DIR
        epr_strip.METRICS_TABLE = PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv"
        epr_strip.main()
        generate_strip_forest_v3()
        return
    if args.zoom_only:
        import generate_publication_figures_v4_filtered_tfs as v4
        import seaborn as sns

        if fig5_paths.SHARED_UNIVERSE:
            shared_universe_axes.apply()

        # In a full run the zooms follow the ATF2 heatmap, whose theme is still active; match it.
        sns.set_theme(style="white", font_scale=0.9)
        generate_k562_macro_zoom(v4)
        return

    bench_root = fig5_paths.BENCH_ROOT / "benchmark_outputs"
    if bench_root.exists():
        shutil.rmtree(bench_root)
    PUBLICATION_TABLE_DIR.mkdir(parents=True, exist_ok=True)

    import benchmark_grn_systems as b
    import compare_regulator_prioritization as compare
    import generate_publication_figures_v4_filtered_tfs as v4
    import generate_system_broad_filtered_publication_figures_v2 as pub_v2
    import plot_publication_v4_epr_strip_forest as epr_strip
    import plot_regulator_prioritization_v2 as reg_v2

    epr_strip.OUT_DIR = V4_FIG_DIR / "5b"
    epr_strip.TABLE_DIR = PUBLICATION_TABLE_DIR
    epr_strip.METRICS_TABLE = PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv"
    install_chronopack_loader(compare)
    if fig5_paths.SHARED_UNIVERSE:
        shared_universe_axes.apply()

    print("== ranked edges (method-specific reference-TF universe)")
    generate_regulator_prioritization_base_tables(compare)
    print("== regulator prioritization v2 tables (top-k shares, rank enrichment)")
    generate_regulator_prioritization_v2_tables(reg_v2)
    print("== per-TF PR curves (5 PR methods)")
    generate_publication_v2_base_tables(pub_v2, b)
    print(f"== per-TF metrics (all {len(b.METHOD_ORDER)} methods; method set {fig5_paths.METHOD_SET})")
    generate_publication_v2_all_method_tables(pub_v2, b)
    print("== v4 figures 5b-5e")
    v4.main()
    print("== v4 EPR strip-forest variants")
    epr_strip.main()
    print("== EPR strip-forests v3 (label placement)")
    generate_strip_forest_v3()
    print("== ATF2 heatmap variant")
    generate_atf2_heatmap_variant()
    print("== K562 ChIP macro-PR zoom (FOS/ELK1-filtered)")
    generate_k562_macro_zoom(v4)
    print(f"Wrote figures to {V4_FIG_DIR}")
    print(f"Wrote tables to {TABLE_DIR}")


if __name__ == "__main__":
    main()
