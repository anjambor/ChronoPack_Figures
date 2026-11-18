#!/usr/bin/env python
"""Side analysis: reference-target enrichment of the gene sets each method ranks, and a toy model of random completion.

1. Per 5b TF (TRRUST: RELA, NFKB1, JUN, CREB1, ATF1, ATF2; K562 ChIP: the same without NFKB1), the rate of reference
   targets in each gene set divided by the rate in all 19,338 evaluation genes (1 = no enrichment):
   ChronoPack_K562's ranked targets (10,202), the BEELINE P < 0.05 TFs + 2000 genes (2,260; all ranked by ChronoPack),
   ChronoPack's top 2,260 targets (same size; its target order is shared by every TF), ChronoPack's targets outside
   the BEELINE set, and the genes ChronoPack does not rank.
2. Toy model of 5b's random-expectation completion (benchmark_grn_systems.evaluate_one, EPR fix): 20,000 genes, 1,000
   targets; method A ranks genes 1-10,000, method B a representative 2,000-gene subset of them. A's genes hold 500
   (not enriched) or 800 (1.6x) targets; inside each method's genes the order is random or better than random
   (score = U(0, 1) + 0.3 for targets). 20 repeats, numpy seed 0.
Writes output/gene_set_enrichment/{enrichment_by_tf.csv, enrichment_median.csv, toy_completion.csv}.
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

OUT = fig5_paths.OUTPUT_DIR / "gene_set_enrichment"
TFS = {"TRRUST": ["RELA", "NFKB1", "JUN", "CREB1", "ATF1", "ATF2"], "K562_ChIP": ["RELA", "JUN", "CREB1", "ATF1", "ATF2"]}


def enrichment() -> pd.DataFrame:
    valid = b.load_valid_targets()
    gts = pub_v2.load_reference_ground_truths(valid)
    order = pd.read_csv(fig5_paths.CHRONOPACK_OUTPUT_DIR / "K562" / "grn.target_scores.csv").target.astype(str)
    cp = [g for g in order if g in valid]
    bee = set(pd.read_csv(fig5_paths.BEELINE_RESULTS_DIR / "GENIE3" / "rankedEdges.csv", sep="\t", usecols=[1]).iloc[:, 0]) & valid
    assert bee <= set(cp), "BEELINE genes not all ranked by ChronoPack"
    sets = {"ChronoPack_all": set(cp), "BEELINE_p05_TFs2000": bee, "ChronoPack_top_same_size": set(cp[:len(bee)]),
            "ChronoPack_not_BEELINE": set(cp) - bee, "not_ranked_by_ChronoPack": valid - set(cp)}
    print({k: len(v) for k, v in sets.items()}, "| BEELINE genes in ChronoPack's top", len(bee), ":",
          len(bee & sets["ChronoPack_top_same_size"]))
    rows = []
    for bench, tfs in TFS.items():
        gt = gts[bench]
        for tf in tfs:
            pos = set(gt.loc[gt.TF == tf, "target"]) & valid
            base = len(pos) / len(valid)
            row = {"benchmark": bench, "TF": tf, "n_targets": len(pos)}
            for name, s in sets.items():
                row[f"{name}_enrichment"] = len(pos & s) / len(s) / base
                row[f"{name}_share_of_targets"] = len(pos & s) / len(pos)
            rows.append(row)
    return pd.DataFrame(rows)


def toy() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    n, k = 20000, 1000
    genes = np.array([f"g{i}" for i in range(n)])

    def score(idx, y):
        e = pd.DataFrame({"TF": "T", "target": genes[idx], "EdgeWeight": np.arange(len(idx), 0, -1, dtype=float)})
        r, _ = b.evaluate_one(e, pd.DataFrame({"TF": "T", "target": genes[y]}), {"T"}, set(genes), curve_points=2)
        return r["AUPRC_ratio"]

    def order(idx, y, boost):
        s = rng.random(len(idx)) + boost * y[idx]
        return idx[np.argsort(-s, kind="mergesort")]

    rows = []
    for label, n_a in [("A not enriched (5% targets)", 500), ("A enriched 1.6x (8% targets)", 800)]:
        for boost in (0.0, 0.3):
            res = {"A": [], "B": []}
            for _ in range(20):
                y = np.zeros(n, bool)
                y[rng.choice(10000, n_a, replace=False)] = True
                y[10000 + rng.choice(10000, k - n_a, replace=False)] = True
                res["A"].append(score(order(np.arange(10000), y, boost), y))
                res["B"].append(score(order(np.arange(2000), y, boost), y))
            rows.append({"A_gene_set": label, "within_set_order": "random" if boost == 0 else "better than random",
                         "AUPRCr_A_10000_genes": np.mean(res["A"]), "AUPRCr_B_2000_subset": np.mean(res["B"])})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    e = enrichment()
    e.to_csv(OUT / "enrichment_by_tf.csv", index=False)
    med = e.groupby("benchmark")[[c for c in e.columns if c.endswith(("_enrichment", "_share_of_targets"))]].median()
    med.to_csv(OUT / "enrichment_median.csv")
    t = toy()
    t.to_csv(OUT / "toy_completion.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(e.round(2).to_string(index=False))
        print(med.round(2).T.to_string())
        print(t.round(2).to_string(index=False))
    print("GENE_SET_ENRICHMENT_OK")


if __name__ == "__main__":
    main()
