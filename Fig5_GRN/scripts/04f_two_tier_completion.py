#!/usr/bin/env python
"""Side analysis: does ChronoPack's 5b lead survive if every comparator gets ChronoPack's candidate gene set?

5b completes each method's ranking with one random tail over all unranked genes. Here each comparator's ranking is
completed in two tiers instead:
  1. the method's own ranked targets, in its own order (unchanged);
  2. the rest of ChronoPack_K562's ranked targets for that TF (ChronoPack's ~10,200 candidate genes), in random order;
  3. every other gene of the 19,338-gene universe, as evaluate_one's random-expectation tail.
So a comparator gets the same gene-set advantage ChronoPack has (predicting for its candidate set before the
background) while keeping its own order inside its ranked genes. Any remaining ChronoPack lead is ordering.

Tier 2 is drawn N_SHUFFLES times (seeds 0..N-1); metrics are averaged. Same edges (the method-specific ranked-edge
table of FIG5_METHOD_SET), references, 5b TF sets, universe and evaluate_one (EPR fix) as 5b. K562 methods only
(ChronoPack_K562's candidate set). ChronoPack_K562 itself is listed unchanged (its candidate set is its ranking).
Writes <BENCH_ROOT>/two_tier_completion/{per_tf_metrics.csv, median_by_method.csv}.
"""
from __future__ import annotations

import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
for p in (SCRIPTS, SCRIPTS / "capsule_lib" / "01_generate_predictions_and_metrics", SCRIPTS / "capsule_lib" / "02_visualize_publication_figures"):
    sys.path.insert(0, str(p))

import fig5_paths  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import benchmark_grn_systems as b  # noqa: E402
import generate_system_broad_filtered_publication_figures_v2 as pub_v2  # noqa: E402
import plot_publication_v4_epr_strip_forest as epr_strip  # noqa: E402

OUT = fig5_paths.BENCH_ROOT / "two_tier_completion"
N_SHUFFLES = 20
CP = "ChronoPack_K562"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    k562 = [m for m in b.METHOD_ORDER if m not in b.THP1_METHODS]
    pub_v2.ALL_METHOD_ORDER = b.METHOD_ORDER
    valid = b.load_valid_targets()
    gts = pub_v2.load_reference_ground_truths(valid)
    requested = {tf for bench in gts for m in k562 for tf in epr_strip.requested_tfs(bench, m)}
    edges = pub_v2.load_ranked_edges_for_v2(valid, requested)
    rows = []
    for bench, gt in gts.items():
        for method in k562:
            for tf in epr_strip.requested_tfs(bench, method):
                e = edges[method]
                e = e[(e["TF"] == tf) & e["target"].isin(valid) & (e["target"] != tf)][["TF", "target", "EdgeWeight"]]
                if e.empty or tf not in set(gt["TF"]):
                    continue
                base, _ = b.evaluate_one(e, gt, {tf}, valid, curve_points=2)
                cp_targets = set(edges[CP].loc[edges[CP]["TF"] == tf, "target"]) & valid
                tier2 = sorted(cp_targets - set(e["target"]) - {tf})
                runs = []
                if method != CP and tier2:
                    floor = float(e["EdgeWeight"].min())
                    for seed in range(N_SHUFFLES):
                        w = floor - 1.0 - np.random.default_rng(seed).permutation(len(tier2)) / len(tier2)
                        e2 = pd.concat([e, pd.DataFrame({"TF": tf, "target": tier2, "EdgeWeight": w})], ignore_index=True)
                        r, _ = b.evaluate_one(e2, gt, {tf}, valid, curve_points=2)
                        runs.append((r["EPR"], r["AUPRC_ratio"]))
                epr2, auprc2 = np.mean(runs, axis=0) if runs else (base["EPR"], base["AUPRC_ratio"])
                rows.append({"benchmark": bench, "method": method, "TF": tf, "k_pos": base["k_pos"],
                             "n_ranked": base["n_pred"], "n_tier2": len(tier2) if method != CP else 0,
                             "EPR_5b": base["EPR"], "AUPRCr_5b": base["AUPRC_ratio"],
                             "EPR_two_tier": epr2, "AUPRCr_two_tier": auprc2})
                print(f"{bench:9s} {method:17s} {tf:6s} AUPRCr {base['AUPRC_ratio']:.2f} -> {auprc2:.2f}", flush=True)
    d = pd.DataFrame(rows)
    d["EPR_5b"] = d["EPR_5b"].clip(upper=pub_v2.EPR_CAP)
    d["EPR_two_tier"] = d["EPR_two_tier"].clip(upper=pub_v2.EPR_CAP)
    d.to_csv(OUT / "per_tf_metrics.csv", index=False)
    med = (d.groupby(["benchmark", "method"], sort=False)
           .agg(n_TFs=("TF", "size"), EPR_5b=("EPR_5b", "median"), EPR_two_tier=("EPR_two_tier", "median"),
                AUPRCr_5b=("AUPRCr_5b", "median"), AUPRCr_two_tier=("AUPRCr_two_tier", "median")).reset_index())
    med.to_csv(OUT / "median_by_method.csv", index=False)
    with pd.option_context("display.width", 200):
        print(med.round(2).to_string(index=False))
    print("TWO_TIER_OK")


if __name__ == "__main__":
    main()
