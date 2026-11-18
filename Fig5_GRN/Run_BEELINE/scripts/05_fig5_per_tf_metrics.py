#!/usr/bin/env python
"""Fig5 5b per-TF EPR / AUPRCr for the r9 BEELINE reruns.

Uses Fig5's own code, read-only (Fig5_GRN's copy, with the EPR top-k fix; fig5_code.py), as
02_reproduce_figures.py does: count_top_predicted_tfs.load_beeline_method
(PPCOR = |corVal| from outFile.txt, the others rankedEdges.csv) restricted to load_reference_tfs, the v2 ground truths (TRRUST, K562 ChIP), the v2 requested TFs
per benchmark and evaluate_one over every valid protein-coding target (targets a method does not rank are completed at
random expectation), EPR capped at 25. Only BEELINE_RESULTS_DIR changes between runs. (Until 2026-11-16 the
original metacell run was also scored and checked against Fig5 v1's table; dropped so nothing here depends on old/.)

Also reports an "own universe" variant (targets restricted to the genes each run ranked), which separates ranking
quality from gene coverage. Fig5 itself uses only the full universe.

Writes output/fig5_metrics/per_tf_metrics.csv and per_tf_comparison.csv.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
os.environ["BEELINE_CONTEXT"] = "K562"  # K562 runs vs Fig5's K562 benchmark
import common as c  # noqa: E402

import fig5_code  # noqa: E402  (Fig5_GRN's scoring code first on sys.path)

OUT_DIR = c.OUTPUT_DIR / "fig5_metrics"
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".mplconfig"))

import benchmark_grn_systems as b  # noqa: E402
import count_top_predicted_tfs as ctp  # noqa: E402
import generate_system_broad_filtered_publication_figures_v2 as pub_v2  # noqa: E402

METHODS = ["GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"]
RUNS = {f"r9_TFs{n}": c.BEELINE_ROOT / "outputs" / c.DATASET_DIR / c.dataset_name(n) for n in c.N_GENES}


def evaluate(edges: pd.DataFrame, gts: dict, universe: set[str]) -> list[dict]:
    rows = []
    for bench, gt in gts.items():
        for tf in pub_v2.method_requested_tfs(bench, "GENIE3"):
            if tf not in set(gt["TF"]) or tf not in set(edges["TF"]):
                continue
            row, _ = b.evaluate_one(edges, gt, {tf}, universe, curve_points=500)
            row["EPR_raw"] = row["EPR"]
            row["EPR"] = min(float(row["EPR"]), pub_v2.EPR_CAP) if pd.notna(row["EPR"]) else row["EPR"]
            rows.append({"benchmark": bench, "TF": tf, **row})
    return rows


def main() -> None:
    valid = b.load_valid_targets()
    gts = pub_v2.load_reference_ground_truths(valid)
    reference_tfs = ctp.load_reference_tfs(valid)
    rows = []
    for run, path in RUNS.items():
        for method in METHODS:
            f = path / method / ("outFile.txt" if method == "PPCOR" else "rankedEdges.csv")
            if not f.exists() or f.stat().st_size == 0:
                print(f"skip {run} {method}: no {f.name}", flush=True)
                continue
            ctp.BEELINE_RESULTS_DIR = path
            edges = ctp.load_beeline_method(method, valid)
            edges = edges[edges["TF"].isin(reference_tfs)]
            edges = edges.sort_values("EdgeWeight", ascending=False, kind="mergesort").reset_index(drop=True)
            own = set(edges["target"])
            for universe_name, universe in [("all_valid_targets", valid), ("own_ranked_targets", own)]:
                for r in evaluate(edges, gts, universe):
                    rows.append({"run": run, "method": method, "universe": universe_name,
                                 "n_ranked_targets": len(own), **r})
            print(f"{run} {method}: {len(edges):,} reference-TF edges, {len(own):,} targets", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(OUT_DIR / "per_tf_metrics.csv", index=False)

    full = res[res.universe == "all_valid_targets"]
    wide = full.pivot_table(index=["benchmark", "TF", "method"], columns="run",
                            values=["EPR", "AUPRC_ratio", "positive_coverage"]).sort_index()
    v4 = pd.read_csv(fig5_code.V4_PER_TF)  # Fig5's ChronoPack (fixed EPR)
    ck = v4[v4.method == "ChronoPack_K562"].set_index(["benchmark", "TF"])[["EPR", "AUPRC_ratio"]]
    wide.to_csv(OUT_DIR / "per_tf_comparison.csv")
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(wide.round(2).to_string())
        print("\nChronoPack_K562 (Fig5):\n" + ck.round(2).to_string())
        own = res[res.universe == "own_ranked_targets"].pivot_table(
            index=["benchmark", "TF", "method"], columns="run", values=["EPR", "AUPRC_ratio"]).sort_index()
        print("\nOwn ranked targets only:\n" + own.round(2).to_string())
        print("\nRanked targets per run:\n" + res.groupby(["run", "method"]).n_ranked_targets.first().unstack().to_string())
    print("FIG5_METRICS_OK")


if __name__ == "__main__":
    main()
