#!/usr/bin/env python
"""Fig5_GRN 04e gene-universe sensitivity, repeated with the r9 BEELINE reruns as the BEELINE methods.

04e (Fig5_GRN/scripts/04e_gene_universe_sensitivity.py) re-evaluates the 5b per-TF EPR / AUPRCr on matched target
universes to separate ranking quality from gene coverage; it was run on the original metacell BEELINE predictions.
This script imports 04e's functions unchanged (Fig5_GRN's copy, with the EPR top-k fix; fig5_code.py) (same edges
for ChronoPack / Replogle, references, 5b TF sets, evaluate_one, oracle ceiling and coverage null) and swaps in the
BEELINE predictions of each run:

- r9_TFs1000, r9_TFs2000: Run_BEELINE outputs, loaded as in 05_fig5_per_tf_metrics.py (Fig5's
  count_top_predicted_tfs loader, reference TFs).
(Until 2026-11-16 the metacell run of Fig5 v1 was also scored and checked against v1's 04e table; dropped so nothing
here depends on old/.)

Universes per run (as in 04e; K562 methods only, so 04e's shared_THP1 is dropped):
- full:        the Fig. 5 universe (19,338 protein-coding targets, random tail for unranked genes).
- hvg:         targets restricted to that run's BEELINE expression genes (GENIE3 / PPCOR targets).
- shared_K562: hvg genes also ranked by ChronoPack_K562.
- full_cp_hvg: full universe, ChronoPack's edges cut to the hvg genes.

Writes output/gene_universe/{per_tf_metrics_by_universe,median_by_universe,positive_coverage_per_tf}.csv.
"""
from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
os.environ["BEELINE_CONTEXT"] = "K562"  # K562 runs vs Fig5's K562 benchmark
import common as c  # noqa: E402

import fig5_code  # noqa: E402  (Fig5_GRN's scoring code first on sys.path)

OUT = c.OUTPUT_DIR / "gene_universe"
OUT.mkdir(parents=True, exist_ok=True)

ge = importlib.import_module("04e_gene_universe_sensitivity")  # sets MPLCONFIGDIR via fig5_paths on import
import count_top_predicted_tfs as ctp  # noqa: E402

b = ge.b
BEELINE = ["GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"]
K562_METHODS = ["ChronoPack_K562", *BEELINE, "Replogle2022_K562"]
RUNS = {f"r9_TFs{n}": c.BEELINE_ROOT / "outputs" / c.DATASET_DIR / c.dataset_name(n) for n in c.N_GENES}


def load_r9_beeline(path: Path, valid: set[str], requested: set[str]) -> dict[str, pd.DataFrame]:
    out = {}
    for method in BEELINE:
        f = path / method / ("outFile.txt" if method == "PPCOR" else "rankedEdges.csv")
        if not f.exists() or f.stat().st_size == 0:
            print(f"  skip {method}: no {f.name}", flush=True)
            continue
        ctp.BEELINE_RESULTS_DIR = path
        e = ctp.load_beeline_method(method, valid)
        e = e[e["TF"].isin(requested)]
        out[method] = (e.sort_values("EdgeWeight", ascending=False, kind="mergesort")
                        .drop_duplicates(["TF", "target"], keep="first")[["TF", "target", "EdgeWeight"]]
                        .reset_index(drop=True))
    return out


def main() -> None:
    b.METHOD_ORDER = K562_METHODS  # ge.pairs / evaluate iterate b.METHOD_ORDER
    valid, gts, base_edges = ge.load_inputs()  # Fig5_GRN tables: ChronoPack, Replogle, Yao
    requested = {tf for bench in gts for m in K562_METHODS for tf in ge.pub_v2.method_requested_tfs(bench, m)}
    reference_tfs = ctp.load_reference_tfs(valid)
    metrics, cov_rows = [], []
    for run, path in RUNS.items():
        print(f"== {run}", flush=True)
        edges = {m: base_edges[m] for m in ["ChronoPack_K562", "Replogle2022_K562"]}
        edges.update(load_r9_beeline(path, valid, requested & reference_tfs))
        b.METHOD_ORDER = [m for m in K562_METHODS if m in edges]
        hvg_src = [m for m in ["GENIE3", "PPCOR"] if m in edges]
        hvg = set().union(*(set(edges[m]["target"]) for m in hvg_src)) & valid
        universes = {"full": valid, "hvg": hvg, "shared_K562": hvg & set(edges["ChronoPack_K562"]["target"]),
                     "full_cp_hvg": valid}
        for k, v in universes.items():
            print(f"  universe {k}: {len(v)} targets", flush=True)
        m = ge.evaluate(edges, gts, universes, hvg)
        m.insert(0, "run", run)
        metrics.append(m)
        for bench, gt in gts.items():
            for method in b.METHOD_ORDER:
                for tf in ge.epr_strip.requested_tfs(bench, method):
                    pos = set(gt.loc[gt["TF"] == tf, "target"]) & valid
                    ranked = set(edges[method].loc[edges[method]["TF"] == tf, "target"])
                    if pos and ranked:
                        cov_rows.append(dict(run=run, benchmark=bench, method=method, TF=tf, n_ranked=len(ranked),
                                             n_pos=len(pos), pos_in_ranked=len(pos & ranked), pos_in_hvg=len(pos & hvg),
                                             n_hvg=len(hvg), n_shared_K562=len(universes["shared_K562"]),
                                             pos_in_shared_K562=len(pos & universes["shared_K562"])))
    m = pd.concat(metrics, ignore_index=True)
    m["EPR_capped"] = m["EPR"].clip(upper=ge.pub_v2.EPR_CAP)
    m.to_csv(OUT / "per_tf_metrics_by_universe.csv", index=False)
    pd.DataFrame(cov_rows).to_csv(OUT / "positive_coverage_per_tf.csv", index=False)

    med = (m[m["k_pos"] > 0].groupby(["benchmark", "run", "universe", "method"])
           .agg(n_TFs=("TF", "size"), median_coverage=("positive_coverage", "median"),
                median_EPR=("EPR_capped", "median"), median_AUPRCr=("AUPRC_ratio", "median"),
                median_null_EPR=("null_EPR", "median"), median_null_AUPRCr=("null_AUPRCr", "median"),
                median_oracle_AUPRCr=("oracle_AUPRCr", "median"))
           .reset_index())
    med.to_csv(OUT / "median_by_universe.csv", index=False)

    with pd.option_context("display.width", 250, "display.max_rows", 500):
        for bench in ["TRRUST", "K562_ChIP"]:
            for metric in ["median_EPR", "median_AUPRCr"]:
                tab = (med[(med.benchmark == bench) & med.universe.isin(["full", "shared_K562"])]
                       .pivot_table(index="method", columns=["universe", "run"], values=metric)
                       .reindex(K562_METHODS))
                print(f"\n== {bench} {metric}\n{tab.round(2).to_string()}")
        sel = m[m.TF.isin(["RELA", "NFKB1", "JUN"]) & m.universe.isin(["full", "shared_K562"])]
        for bench in ["TRRUST", "K562_ChIP"]:
            for metric in ["EPR_capped", "AUPRC_ratio"]:
                tab = (sel[sel.benchmark == bench]
                       .pivot_table(index=["TF", "method"], columns=["universe", "run"], values=metric))
                print(f"\n== {bench} per-TF {metric}\n{tab.round(2).to_string()}")
        print("\nShared_K562 universe size and positives:\n" + pd.DataFrame(cov_rows)
              .query("method == 'ChronoPack_K562'").pivot_table(index=["benchmark", "TF"], columns="run",
                                                                values=["n_shared_K562", "pos_in_shared_K562"]).to_string())
    print("GENE_UNIVERSE_OK")


if __name__ == "__main__":
    main()
