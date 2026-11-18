#!/usr/bin/env python
"""Sensitivity of ChronoPack's global GRN ranking to how the TF prior and target score are combined
(side analysis; Fig5 outputs and the ChronoPack implementation are not changed).

Inputs are unchanged: r_TF = priority_act share from the r9 OmniPath TF ranking (tf_ranking.csv) and
g = 1 - lfdr_bum of the tradeSeq association p-values (target_scores.scores.csv). Only the combination differs:

  original     r_TF * g
  percentile   pct(r_TF) * pct(g)            percentile ranks in (0, 1], average rank for ties
  zsum         z(r_TF) + z(g)                z-scores of the raw values
  sigmoid      sigmoid(z(r_TF)) * sigmoid(z(g))  (fixed logistic, no tuned parameters)

TF statistics are taken over the benchmark's universe TFs (TRRUST 45, K562 ChIP 28: TFs predicted by both
ChronoPack contexts and GENIE3/GRNBOOST2/PPCOR/SINCERITIES), gene statistics over ChronoPack's scored genes.
Ties: within a TF, equal combined scores keep ChronoPack's own secondary key (the tradeSeq Wald statistic);
remaining exact ties are taken in random order (expected values), as in 04.

Metrics, all over the full universe TFs x valid targets (variant B):
  pooled EPR / AUPRCr   all TFs' edges in one ranking, all reference edges of the universe TFs (b.evaluate_one;
                        expected TP at k for EPR; AP from the ranking ordered by score, Wald, TF prior)
  TF-balanced recall    the 04 metric for the Fig5 signaling TFs, with and without RELA/NFKB1 rewarded
  edge distribution     expected edges per TF among the top 1,000 / 10,000 edges
GENIE3 and GRNBOOST2 are included for reference. Writes output/score_combination/.
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
for path in (SCRIPTS_DIR, SCRIPTS_DIR / "capsule_lib" / "01_generate_predictions_and_metrics", SCRIPTS_DIR / "capsule_lib" / "02_visualize_publication_figures"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import fig5_paths

os.environ.setdefault("MPLCONFIGDIR", str(fig5_paths.OUTPUT_DIR / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import benchmark_grn_systems as b
import compare_regulator_prioritization as compare
import generate_system_broad_filtered_publication_figures_v2 as pub_v2

spec = importlib.util.spec_from_file_location("g04", SCRIPTS_DIR / "04_tf_balanced_global_recall.py")
g04 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g04)

OUT = fig5_paths.OUTPUT_DIR / "score_combination"
SCHEMES = {"original": "r_TF x (1-lfdr)", "percentile": "percentile product", "zsum": "z-score sum", "sigmoid": "sigmoid(z) product"}
REFERENCE = ["GENIE3", "GRNBOOST2"]
COLORS = {"original": "#111111", "percentile": "#d95f02", "zsum": "#1b9e77", "sigmoid": "#7570b3", "GENIE3": "#8c8c8c", "GRNBOOST2": "#bdbdbd", "Random": "#cccccc"}


def pct(x: pd.Series) -> pd.Series:
    return x.rank(method="average", ascending=True) / len(x)


def z(x: pd.Series) -> pd.Series:
    return (x - x.mean()) / x.std(ddof=0)


def combined_edges(prior: pd.Series, genes: pd.DataFrame, tfs: list[str], scheme: str) -> pd.DataFrame:
    r = prior.reindex(tfs).astype(float)
    g = genes.set_index("gene")["score"].astype(float)
    if scheme == "original":
        a, c, op = r, g, "mul"
    elif scheme == "percentile":
        a, c, op = pct(r), pct(g), "mul"
    elif scheme == "zsum":
        a, c, op = z(r), z(g), "add"
    elif scheme == "sigmoid":
        a, c, op = 1 / (1 + np.exp(-z(r))), 1 / (1 + np.exp(-z(g))), "mul"
    else:
        raise ValueError(scheme)
    av, cv = a.to_numpy(), c.reindex(genes["gene"]).to_numpy()
    key1 = (np.add.outer(av, cv) if op == "add" else np.outer(av, cv)).ravel()
    df = pd.DataFrame({"TF": np.repeat(tfs, len(cv)), "target": np.tile(genes["gene"].to_numpy(), len(tfs)),
                       "key1": key1, "key2": np.tile(genes["waldStat_1"].to_numpy(float), len(tfs)),
                       "prior": np.repeat(r.to_numpy(), len(cv))})
    return df[df["TF"] != df["target"]].reset_index(drop=True)


def pooled_metrics(edges: pd.DataFrame, gt: pd.DataFrame, tfs: list[str], valid: set[str], n_universe: int,
                   positives: set, pos_by_tf: dict) -> tuple[dict, pd.DataFrame]:
    """Pooled EPR (expected, ties random) and AUPRCr (b.evaluate_one on the score / Wald / prior order)."""
    e = edges.copy()
    if "prior" not in e:
        e["prior"] = 0.0
    order = np.lexsort((-e["prior"].to_numpy(), -e["key2"].fillna(-np.inf).to_numpy(), -e["key1"].to_numpy()))
    e = e.iloc[order].reset_index(drop=True)
    e["EdgeWeight"] = np.arange(len(e), 0, -1, dtype=float)
    row, curve = b.evaluate_one(e[["TF", "target", "EdgeWeight"]], gt, set(tfs), valid, curve_points=400)
    k = int(row["k_pos"])
    rk = g04.Ranking(edges, n_universe)
    tp = rk.recall(positives, pos_by_tf, tfs, np.array([float(k)])).mul(pd.Series(pos_by_tf).reindex(tfs), axis=0).sum().iloc[0]
    epr = (tp / k) / row["rho"]
    return {"k_pos": k, "EPR": float(epr), "EPR_evaluate_one": float(row["EPR"]), "AUPRCr": float(row["AUPRC_ratio"])}, curve


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(g04.SEED)
    valid = b.load_valid_targets()
    gts = pub_v2.load_reference_ground_truths(valid)
    gt_tfs = set().union(*[set(g["TF"]) for g in gts.values()])
    coexp = {}
    for m in g04.COEXP:
        e = compare.c.load_beeline_method(m, valid)
        coexp[m] = e[e["TF"].isin(gt_tfs)].rename(columns={"EdgeWeight": "key1"}).assign(key2=0.0)[["TF", "target", "key1", "key2"]]
    comps = {ctx: g04.chronopack_components(ctx, valid)[:2] for ctx in g04.CONTEXT_TFS}

    rows, dist_rows, pr_rows, curve_rows = [], [], [], []
    ks = g04.K_GRID.astype(float)
    for bench in g04.BENCHMARKS:
        gt = gts[bench]
        avail = [set(e["TF"]) for e in coexp.values()] + [set(c[0].index) for c in comps.values()]
        tfs = sorted(set(gt["TF"]).intersection(*avail))
        n_universe = len(tfs) * len(valid)
        gt_u = gt[gt["TF"].isin(tfs) & gt["target"].isin(valid) & (gt["TF"] != gt["target"])].drop_duplicates(["TF", "target"])
        positives = set(zip(gt_u["TF"], gt_u["target"]))
        pos_by_tf = {t: int(gt_u["TF"].eq(t).sum()) for t in tfs}
        for ctx, (prior, genes) in comps.items():
            signaling = [t for t in g04.CONTEXT_TFS[ctx] if t in tfs and pos_by_tf[t] >= g04.MIN_TARGETS]
            sets = {"all_signaling": signaling, "without_RELA_NFKB1": [t for t in signaling if t not in g04.EXCLUDED_DOMINANT]}
            rankings = {s: combined_edges(prior, genes, tfs, s) for s in SCHEMES}
            rankings.update({m: coexp[m][coexp[m]["TF"].isin(tfs)] for m in REFERENCE})
            for name, e in rankings.items():
                e = e[e["target"].isin(valid)]
                pooled, curve = pooled_metrics(e, gt, tfs, valid, n_universe, positives, pos_by_tf)
                rk = g04.Ranking(e, n_universe)
                rec = rk.recall(positives, pos_by_tf, tfs, ks)
                share = rk.edge_share(tfs, len(valid), np.array([1000.0, 10000.0]))
                row = {"benchmark": bench, "context": ctx, "scheme": name, "n_TFs_universe": len(tfs), **pooled}
                for set_name, tfl in sets.items():
                    c = rec.loc[tfl].mean(axis=0)
                    sm = g04.summarise(bench, ctx, "all", set_name, name, tfl, rec, c, n_universe, rng)
                    row.update({f"{set_name}_logAUC_vs_random": sm["logAUC_vs_random"], f"{set_name}_R@1000": sm["R@1000"],
                                f"{set_name}_R@10000": sm["R@10000"], f"{set_name}_nTFs_recall>=0.05@10000": sm[f"n_TFs_recall>={g04.MEANINGFUL_RECALL:g}@10000"]})
                    for kk, v in c.items():
                        curve_rows.append({"benchmark": bench, "context": ctx, "scheme": name, "rewarded": set_name, "K": kk, "R_bal": v})
                for kcol, K in ((1000.0, 1000), (10000.0, 10000)):
                    s = share[kcol] / K
                    row[f"top{K}_n_TFs_with_edges"] = int((share[kcol] >= 1).sum())
                    row[f"top{K}_effective_n_TFs"] = float(1 / (s ** 2).sum())
                    row[f"top{K}_RELA_NFKB1_share"] = float(s.reindex(["RELA", "NFKB1"]).fillna(0).sum())
                    row[f"top{K}_signaling_share"] = float(s.reindex(signaling).sum())
                    top3 = s.sort_values(ascending=False).head(3)
                    row[f"top{K}_largest_TFs"] = ", ".join(f"{t} {v:.0%}" for t, v in top3.items())
                    for tf, v in s.items():
                        dist_rows.append({"benchmark": bench, "context": ctx, "scheme": name, "budget": K, "TF": tf, "share": v, "signaling": tf in signaling})
                rows.append(row)
                c = curve.assign(benchmark=bench, context=ctx, scheme=name)
                pr_rows.append(c)
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "summary.csv", index=False)
    pd.DataFrame(dist_rows).to_csv(OUT / "edge_distribution_by_tf.csv", index=False)
    pd.concat(pr_rows).to_csv(OUT / "pooled_pr_curves.csv", index=False)
    curves = pd.DataFrame(curve_rows)
    curves.to_csv(OUT / "tf_balanced_curves.csv", index=False)

    pd.set_option("display.width", 260)
    cols = ["scheme", "EPR", "AUPRCr", "all_signaling_logAUC_vs_random", "without_RELA_NFKB1_logAUC_vs_random",
            "all_signaling_R@10000", "without_RELA_NFKB1_R@10000", "top1000_effective_n_TFs", "top1000_RELA_NFKB1_share",
            "top10000_effective_n_TFs", "top10000_RELA_NFKB1_share", "top10000_largest_TFs"]
    for (bench, ctx), g in res.groupby(["benchmark", "context"], sort=False):
        print(f"\n== {bench} | {ctx} | universe {g.n_TFs_universe.iloc[0]} TFs (EPR check: expected vs evaluate_one max diff {np.abs(g.EPR - g.EPR_evaluate_one).max():.3f}) ==")
        print(g[cols].round(3).to_string(index=False))
    plots(res, curves, pd.concat(pr_rows))
    print(f"\nWrote {OUT}")


def plots(res: pd.DataFrame, curves: pd.DataFrame, pr: pd.DataFrame) -> None:
    plt.rcParams.update({"pdf.fonttype": 42, "font.size": 8})
    combos = list(res.groupby(["benchmark", "context"], sort=False).groups)
    fig, axes = plt.subplots(len(combos), 3, figsize=(15, 3.3 * len(combos)))
    for i, (bench, ctx) in enumerate(combos):
        for j, rew in enumerate(["all_signaling", "without_RELA_NFKB1"]):
            ax = axes[i, j]
            d = curves[(curves.benchmark == bench) & (curves.context == ctx) & (curves.rewarded == rew)]
            for name in list(SCHEMES) + REFERENCE:
                dv = d[d.scheme == name].sort_values("K")
                ax.plot(dv.K, dv.R_bal, color=COLORS[name], lw=1.8 if name in SCHEMES else 1.1, ls="-" if name in SCHEMES else "--",
                        label=SCHEMES.get(name, name))
            ax.set_xscale("log")
            ax.set_title(f"{bench} | {ctx} | TF-balanced recall, {'all signaling TFs' if rew == 'all_signaling' else 'without RELA/NFKB1'}", loc="left", fontsize=7.5)
            ax.grid(True, color="#eeeeee")
            ax.set_xlabel("Global edge budget K")
            if j == 0:
                ax.set_ylabel("TF-balanced recall")
            if i == 0 and j == 0:
                ax.legend(fontsize=7, frameon=False)
        ax = axes[i, 2]
        d = pr[(pr.benchmark == bench) & (pr.context == ctx)]
        for name in list(SCHEMES) + REFERENCE:
            dv = d[(d.scheme == name) & (~d.is_expectation)].sort_values("recall")
            ax.plot(dv.recall, dv.precision, color=COLORS[name], lw=1.8 if name in SCHEMES else 1.1, ls="-" if name in SCHEMES else "--")
        ax.set_xscale("log")
        ax.set_xlabel("Pooled recall (all universe TFs)")
        ax.set_ylabel("Precision")
        ax.set_title(f"{bench} | {ctx} | pooled global PR (ranked edges)", loc="left", fontsize=7.5)
        ax.grid(True, color="#eeeeee")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"score_combination_curves.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
