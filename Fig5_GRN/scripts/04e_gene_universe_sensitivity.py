#!/usr/bin/env python
"""Side analysis: how much of the per-TF EPR / AUPRCr gap in Fig. 5b comes from target-gene coverage?

Fig. 5 scores every method against the same broad universe (all filtered protein-coding genes, 19,338 targets per
TF) and completes each method's ranking with a random tail. The BEELINE methods rank only the ~2.9k genes of their
expression input (GRNBOOST2 even fewer, ~400 per TF), whereas ChronoPack ranks ~10k tradeSeq-fit genes and the
perturbation screens 7.7k-13k. Reference positives outside a method's ranked set land in the random tail.

Re-evaluates the 5b per-TF metrics (same edges, references, TF sets and evaluate_one) under:

- full:        the Fig. 5 universe (reproduces the published table).
- hvg:         targets restricted to the BEELINE expression gene set (GENIE3/PPCOR targets; 2,918 genes in Fig5_GRN,
               2,338 in Fig5_GRN_v2).
- shared_K562: hvg genes also ranked by ChronoPack_K562 (2,478), where ChronoPack_K562 and GENIE3/PPCOR/SINCERITIES
               each rank (nearly) every gene.
- shared_THP1: hvg genes also ranked by ChronoPack_THP1 (1,910).
- full_cp_hvg: full universe, but ChronoPack's edges cut to hvg targets (ChronoPack predicting only for the
               BEELINE genes; everything else becomes random tail, as for the BEELINE methods).
and, in the full universe, two coverage references per method: the oracle ceiling (AUPRCr if every positive inside the
method's ranked set were ranked first) and the coverage null (mean AUPRCr / EPR over N_SHUFFLES random orders of the
method's own ranked set, i.e. what knowing only *which* genes a method ranks is worth).

Nothing in Fig. 5 is changed. Outputs: output/gene_universe/.
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

OUT = fig5_paths.BENCH_ROOT / "gene_universe"  # Fig5_GRN_v3: output/gene_universe or output/with_prior/gene_universe
PUBLISHED = fig5_paths.PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv"
CHRONOPACK = ["ChronoPack_K562", "ChronoPack_THP1"]
N_SHUFFLES = 20


def load_inputs():
    # same set-up as 02_reproduce_figures.generate_publication_v2_all_method_tables
    pub_v2.ALL_METHOD_ORDER = b.METHOD_ORDER
    pub_v2.METHOD_CONTEXT = {m: ("THP1" if m in {"ChronoPack_THP1", "Yao2023_THP1"} else "K562") for m in b.METHOD_ORDER}
    valid = b.load_valid_targets()
    gts = pub_v2.load_reference_ground_truths(valid)
    requested = {tf for bench in gts for m in b.METHOD_ORDER for tf in pub_v2.method_requested_tfs(bench, m)}
    edges = pub_v2.load_ranked_edges_for_v2(valid, requested)
    return valid, gts, edges


def oracle_auprcr(n_pred: int, n_universe: int, k_pos: int, t: int) -> float:
    """AUPRCr when the t positives inside the ranked set are ranked 1..t and the rest are random tail."""
    ap = b.ap_random_completion(np.arange(t), n_pred, n_universe, k_pos)
    return ap / (k_pos / n_universe)


def coverage_null(e: pd.DataFrame, gt: pd.DataFrame, tf: str, universe: set[str]) -> tuple[float, float]:
    rng = np.random.default_rng(0)
    out = []
    for _ in range(N_SHUFFLES):
        s = e.assign(EdgeWeight=rng.permutation(len(e)).astype(float))
        r, _ = b.evaluate_one(s, gt, {tf}, universe, curve_points=2)
        out.append((r["EPR"], r["AUPRC_ratio"]))
    return tuple(np.mean(out, axis=0))


def evaluate(edges, gts, universes, cp_hvg):
    rows = []
    for (bench, method, tf) in pairs(gts, edges):
        for uname, universe in universes.items():
            if uname == "full_cp_hvg":
                if method not in CHRONOPACK:
                    continue
                e = edges[method][edges[method]["target"].isin(cp_hvg)]
            else:
                e = edges[method]
            row, _ = b.evaluate_one(e, gts[bench], {tf}, universe, curve_points=50)
            row.update(benchmark=bench, method=method, TF=tf, universe=uname, n_universe_targets=len(universe))
            if uname == "full" and row["k_pos"]:
                t = int(row["TP_emitted"])
                row["oracle_AUPRCr"] = oracle_auprcr(int(row["n_pred"]), int(row["N"]), int(row["k_pos"]), t)
                row["null_EPR"], row["null_AUPRCr"] = coverage_null(e[e["TF"] == tf], gts[bench], tf, universe)
            rows.append(row)
    return pd.DataFrame(rows)


def pairs(gts, edges):
    for bench, gt in gts.items():
        gt_tfs = set(gt["TF"])
        for method in b.METHOD_ORDER:
            for tf in epr_strip.requested_tfs(bench, method):  # the 5b TF set (FOS/ELK1 dropped as in 5b)
                if tf in gt_tfs and tf in set(edges[method]["TF"]):
                    yield bench, method, tf


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    valid, gts, edges = load_inputs()
    hvg = (set(edges["GENIE3"]["target"]) | set(edges["PPCOR"]["target"])) & valid
    universes = {
        "full": valid,
        "hvg": hvg,
        "shared_K562": hvg & set(edges["ChronoPack_K562"]["target"]),
        "shared_THP1": hvg & set(edges["ChronoPack_THP1"]["target"]),
        "full_cp_hvg": valid,
    }
    for k, v in universes.items():
        print(f"universe {k}: {len(v)} targets")

    # coverage: ranked targets per method and where the reference positives sit
    cov = []
    for bench, gt in gts.items():
        for method in b.METHOD_ORDER:
            for tf in epr_strip.requested_tfs(bench, method):
                pos = set(gt.loc[gt["TF"] == tf, "target"]) & valid
                ranked = set(edges[method].loc[edges[method]["TF"] == tf, "target"])
                if not pos or not ranked:
                    continue
                cov.append(dict(benchmark=bench, method=method, TF=tf, n_ranked=len(ranked), n_pos=len(pos),
                                pos_in_ranked=len(pos & ranked), pos_in_hvg=len(pos & hvg),
                                pos_rate_full=len(pos) / len(valid), pos_rate_hvg=len(pos & hvg) / len(hvg)))
    cov = pd.DataFrame(cov)
    cov.to_csv(OUT / "positive_coverage_per_tf.csv", index=False)

    m = evaluate(edges, gts, universes, hvg)
    m.to_csv(OUT / "per_tf_metrics_by_universe.csv", index=False)

    # the full universe must reproduce the published 5b table (EPR there is capped at 25)
    pub = pd.read_csv(PUBLISHED).set_index(["benchmark", "method", "TF"])
    full = m[m["universe"] == "full"].set_index(["benchmark", "method", "TF"])
    j = full.join(pub[["EPR_raw", "AUPRC_ratio"]], rsuffix="_pub", how="inner")
    assert len(j) == len(full), "some 5b rows missing from the published table"
    assert np.allclose(j["EPR"], j["EPR_raw"], equal_nan=True) and np.allclose(j["AUPRC_ratio"], j["AUPRC_ratio_pub"], equal_nan=True)
    print(f"full universe reproduces the published per-TF table ({len(j)} rows)")

    m["EPR_capped"] = m["EPR"].clip(upper=pub_v2.EPR_CAP)
    med = (m[m["k_pos"] > 0].groupby(["benchmark", "universe", "method"])
           .agg(n_TFs=("TF", "size"), median_k_pos=("k_pos", "median"), median_coverage=("positive_coverage", "median"),
                median_EPR=("EPR_capped", "median"), median_AUPRCr=("AUPRC_ratio", "median"),
                median_oracle_AUPRCr=("oracle_AUPRCr", "median"),
                median_null_EPR=("null_EPR", "median"), median_null_AUPRCr=("null_AUPRCr", "median"))
           .reset_index())
    med.to_csv(OUT / "median_by_universe.csv", index=False)

    pd.set_option("display.width", 220)
    order = {u: i for i, u in enumerate(universes)}
    for bench in ["TRRUST", "K562_ChIP"]:
        sub = med[med["benchmark"] == bench].copy()
        for metric in ["median_EPR", "median_AUPRCr"]:
            tab = sub.pivot(index="method", columns="universe", values=metric).reindex(b.METHOD_ORDER)
            tab = tab[[u for u in universes if u in tab.columns]]
            print(f"\n== {bench} {metric}\n{tab.round(2).to_string()}")
        tab = sub[sub["universe"] == "full"].set_index("method").reindex(b.METHOD_ORDER)
        print(f"\n== {bench} full-universe coverage and oracle ceiling\n"
              f"{tab[['median_coverage', 'median_EPR', 'median_null_EPR', 'median_AUPRCr', 'median_null_AUPRCr', 'median_oracle_AUPRCr']].round(2).to_string()}")
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
