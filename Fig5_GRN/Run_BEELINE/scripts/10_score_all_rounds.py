#!/usr/bin/env python
"""Fig5 5b per-TF metrics and the 04e matched-gene analysis for every BEELINE round, both contexts.

Rounds (TFs + 1000 and TFs + 2000 each; K562 and THP1):
  r9_10pct      BEELINE 10% expression filter, all TFs + top-N variance genes (-p 0)    <ctx>-r9-sc-TFs<N>
  r9_nofilter   no expression filter (ChronoPack's), -p 0                               <ctx>-r9-sc-cpfilt-TFs<N>
  r9_p05        no expression filter, uncorrected GAM P < 0.05, significant TFs + top-N <ctx>-r9-sc-cpfilt-p05-TFs<N>
                variance significant genes (the BEELINE paper's selection, looser cut-off)
(Until 2026-11-16 also K562 metacell_Fig5, the BEELINE run Fig5 v1 used; dropped so nothing here depends on old/.)

Scoring is Fig5's own code, imported read-only, as in steps 05 / 06 (Fig5_GRN's copy, with the EPR top-k fix;
fig5_code.py): 04e's load_inputs (ChronoPack, Replogle and
Yao edges from Fig5's ranked-edge table), count_top_predicted_tfs.load_beeline_method for BEELINE (PPCOR = |corVal|
from outFile.txt, the others rankedEdges.csv), restricted to reference TFs, and 04e's evaluate (evaluate_one with the
5b TF sets, oracle ceiling and coverage null in the full universe). A THP1 BEELINE method is scored like the THP1
methods of Fig5 (ChronoPack_THP1, Yao2023_THP1): THP1 TF set (adds REL and IRF3; FOS dropped), same references.
THP1 PPCOR is not applicable (787 cells < genes; README) and is not scored.

Universes per run (04e definitions; <ctx> = the run's context):
  full          the Fig5 universe (19,338 protein-coding targets, random tail for unranked genes)
  hvg           the run's BEELINE expression genes
  shared        hvg genes also ranked by ChronoPack_<ctx> (04e's shared_K562 / shared_THP1)
  full_cp_hvg   full universe, ChronoPack_<ctx> edges cut to the hvg genes

Checks: the K562 r9_10pct rows must reproduce step 06's table (output/gene_universe; same code and fix);
ChronoPack / Replogle / Yao full-universe rows must equal Fig5_GRN's per-TF table (fixed EPR).
Writes output/rounds/{per_tf_metrics.csv, median_by_universe.csv, runs.csv}.
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

import fig5_code  # noqa: E402  (Fig5_GRN's scoring code first on sys.path)

OUT = c.OUTPUT_DIR / "rounds"
OUT.mkdir(parents=True, exist_ok=True)

ge = importlib.import_module("04e_gene_universe_sensitivity")
import count_top_predicted_tfs as ctp  # noqa: E402

b = ge.b
BEELINE = ["GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"]
FIG5_ORDER = list(b.METHOD_ORDER)
OTHERS = {"K562": ["ChronoPack_K562", "Replogle2022_K562"], "THP1": ["ChronoPack_THP1", "Yao2023_THP1"]}
ROUNDS = {"r9_10pct": "", "r9_nofilter": "-cpfilt", "r9_p05": "-cpfilt-p05"}
STEP06 = c.OUTPUT_DIR / "gene_universe" / "per_tf_metrics_by_universe.csv"
PUBLISHED = fig5_code.V4_PER_TF


def key(method: str, ctx: str) -> str:
    """Method key: K562 BEELINE keeps Fig5's names; THP1 BEELINE gets a _THP1 suffix."""
    return method if ctx == "K562" else f"{method}_THP1"


def set_contexts(keys: list[str]) -> None:
    """Make Fig5's TF-set helpers treat the *_THP1 keys (and ChronoPack_THP1 / Yao2023_THP1) as THP1 methods."""
    thp1 = {k for k in keys if k.endswith("_THP1")}
    ge.epr_strip.THP1_METHODS = set(ge.epr_strip.THP1_METHODS) | thp1
    ge.pub_v2.METHOD_CONTEXT = {**ge.pub_v2.METHOD_CONTEXT, **{k: ("THP1" if k in ge.epr_strip.THP1_METHODS else "K562") for k in keys}}


def load_beeline(path: Path, ctx: str, valid: set[str], tfs: set[str]) -> dict[str, pd.DataFrame]:
    out = {}
    for method in BEELINE:
        if ctx == "THP1" and method == "PPCOR":
            continue  # not applicable (more genes than cells)
        f = path / method / ("outFile.txt" if method == "PPCOR" else "rankedEdges.csv")
        if not f.exists() or f.stat().st_size == 0:
            print(f"  SKIP {path.name} {method}: no {f.name} (run not finished)", flush=True)
            continue
        ctp.BEELINE_RESULTS_DIR = path
        e = ctp.load_beeline_method(method, valid)
        e = e[e["TF"].isin(tfs)]
        out[key(method, ctx)] = (e.sort_values("EdgeWeight", ascending=False, kind="mergesort")
                                 .drop_duplicates(["TF", "target"], keep="first")[["TF", "target", "EdgeWeight"]]
                                 .reset_index(drop=True))
    return out


def main() -> None:
    b.METHOD_ORDER = FIG5_ORDER
    valid, gts, base_edges = ge.load_inputs()
    reference_tfs = ctp.load_reference_tfs(valid)
    all_keys = FIG5_ORDER + [key(m, "THP1") for m in BEELINE if m != "PPCOR"]
    set_contexts(all_keys)
    requested = {tf for bench in gts for m in all_keys for tf in ge.pub_v2.method_requested_tfs(bench, m)}
    tfs = requested & reference_tfs

    runs = [(ctx, rnd, n, c.BEELINE_ROOT / "outputs" / c.DATASET_DIR / f"{ctx}-r9-sc{tag}-TFs{n}")
             for ctx in ("K562", "THP1") for rnd, tag in ROUNDS.items() for n in c.N_GENES]
    metrics, run_rows = [], []
    for ctx, rnd, n, path in runs:
        name = rnd if n is None else f"{rnd}_TFs{n}"
        print(f"== {ctx} {name}", flush=True)
        edges = {m: base_edges[m] for m in OTHERS[ctx]}
        edges.update(load_beeline(path, ctx, valid, tfs))
        genes = pd.read_csv(c.BEELINE_ROOT / "BEELINE-data" / "inputs" / c.DATASET_DIR / path.name / "ExpressionData.csv",
                            index_col=0, usecols=[0]).index
        hvg = set(genes) & valid
        g3 = set(edges[key("GENIE3", ctx)]["target"])
        assert g3 <= hvg and len(hvg - g3) <= 5, f"GENIE3 targets != expression genes ({len(hvg - g3)} missing)"
        cp = f"ChronoPack_{ctx}"
        b.METHOD_ORDER = [m for m in all_keys if m in edges]
        universes = {"full": valid, "hvg": hvg, "shared": hvg & set(edges[cp]["target"]), "full_cp_hvg": valid}
        m = ge.evaluate(edges, gts, universes, hvg)
        m = m[~((m.universe == "full_cp_hvg") & (m.method != cp))]  # 04e: ChronoPack of this context only
        m.insert(0, "run", name)
        m.insert(0, "TFsN", n)
        m.insert(0, "round", rnd)
        m.insert(0, "context", ctx)
        metrics.append(m)
        run_rows.append({"context": ctx, "run": name, "round": rnd, "TFsN": n, "dataset": path.name,
                         "n_hvg": len(hvg), "n_shared": len(universes["shared"]),
                         **{f"n_targets_{k}": int(edges[k]["target"].nunique()) for k in b.METHOD_ORDER}})
        print("  " + ", ".join(f"{k} {len(v)}" for k, v in universes.items()), flush=True)
    m = pd.concat(metrics, ignore_index=True)
    m["EPR_capped"] = m["EPR"].clip(upper=ge.pub_v2.EPR_CAP)
    m.to_csv(OUT / "per_tf_metrics.csv", index=False)
    pd.DataFrame(run_rows).to_csv(OUT / "runs.csv", index=False)

    # checks
    ref = pd.read_csv(STEP06).replace({"universe": {"shared_K562": "shared"},
                                      "run": {f"r9_TFs{n}": f"r9_10pct_TFs{n}" for n in c.N_GENES}})
    ref = ref[ref.universe.isin(["full", "hvg", "shared", "full_cp_hvg"])]
    mine = m[m.context.eq("K562") & m.run.isin(ref.run.unique())]
    j = mine.merge(ref, on=["run", "benchmark", "method", "TF", "universe"], suffixes=("", "_06"))
    assert len(j) == len(ref) == len(mine), (len(j), len(ref), len(mine))
    for col in ["EPR", "AUPRC_ratio"]:
        assert np.allclose(j[col], j[f"{col}_06"], equal_nan=True), col
    print(f"K562 r9_10pct rows reproduce step 06 ({len(j)} rows)")
    pub = pd.read_csv(PUBLISHED).set_index(["benchmark", "method", "TF"])
    oth = m[m.universe.eq("full") & m.method.isin(OTHERS["K562"] + OTHERS["THP1"])].set_index(["benchmark", "method", "TF"])
    jj = oth.join(pub[["EPR_raw", "AUPRC_ratio"]], rsuffix="_pub", how="inner")
    assert len(jj) == len(oth) and np.allclose(jj["EPR"], jj["EPR_raw"]) and np.allclose(jj["AUPRC_ratio"], jj["AUPRC_ratio_pub"])
    print(f"ChronoPack / Replogle / Yao full-universe rows equal Fig5_GRN per-TF table ({len(jj)} rows)")

    med = (m[m["k_pos"] > 0].groupby(["context", "benchmark", "round", "TFsN", "run", "universe", "method"], dropna=False)
           .agg(n_TFs=("TF", "size"), median_coverage=("positive_coverage", "median"),
                median_EPR=("EPR_capped", "median"), median_AUPRCr=("AUPRC_ratio", "median"),
                median_null_EPR=("null_EPR", "median"), median_null_AUPRCr=("null_AUPRCr", "median"),
                median_oracle_AUPRCr=("oracle_AUPRCr", "median"))
           .reset_index())
    med.to_csv(OUT / "median_by_universe.csv", index=False)
    with pd.option_context("display.width", 260, "display.max_rows", 500, "display.max_columns", 40):
        for bench in ["TRRUST", "K562_ChIP"]:
            for metric in ["median_EPR", "median_AUPRCr"]:
                for uni in ["full", "shared"]:
                    tab = (med[(med.benchmark == bench) & (med.universe == uni)]
                           .pivot_table(index=["context", "method"], columns="run", values=metric))
                    print(f"\n== {bench} {metric} {uni}\n{tab.round(2).to_string()}")
    print("ROUNDS_SCORING_OK")


if __name__ == "__main__":
    main()
