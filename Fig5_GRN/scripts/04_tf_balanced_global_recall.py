#!/usr/bin/env python
"""Signaling-TF-balanced global edge recovery (side analysis; Fig5 outputs are not changed).

Question: does a method jointly rank stimulus-relevant TFs and their known targets near the top of
its *global* edge list, or is pooled global enrichment explained by the annotation density of
RELA/NFKB1 in TRRUST?

Universe (per benchmark) -- the common-TF universe of the earlier global analysis: TFs that both
ChronoPack contexts and the four co-expression methods (GENIE3, GRNBOOST2, PPCOR, SINCERITIES)
predict for and that the reference covers (TRRUST: 45, K562 ChIP: 28), x all valid targets
(b.load_valid_targets, as in every Fig5 evaluation). All edges in this universe are ranked; the
predicted network is NOT restricted to signaling TFs.

Rewarded TFs -- the signaling TFs annotated in Fig5 (K562_EVAL_TFS / THP1_EVAL_TFS of
plot_regulator_prioritization_v2: RELA, NFKB1, JUN, ELK1, CREB1, ATF1, FOS, ATF2; THP1 adds REL, IRF3),
fixed before looking at any result here, kept if in the universe with >= MIN_TARGETS reference targets.

Metric -- TF-balanced recall R_bal(K) = mean over rewarded TFs t of (expected known targets of t among
the top K edges) / (known targets of t). Ties: edges with exactly equal scores (after a method's own
documented secondary key, e.g. ChronoPack's tradeSeq Wald statistic) are taken in random order, so each
contributes its expected value; edges a method does not rank at all form a final randomly ordered tail.
Random global ranking: R_t(K) = K / N for every TF (N = universe size). Summaries: R_bal at fixed
budgets and logAUC = mean R_bal over log-spaced K in [10, 10^4] (also as a ratio to random).
Uncertainty: bootstrap over rewarded TFs (few TFs, so intervals are wide and only indicative).

Methods -- co-expression baselines (BEELINE outputs on K562 data, raw scores via the capsule loader;
the same rankings are compared with both ChronoPack contexts). Perturbation baselines (Replogle, Yao)
cover only ~10 TFs and are excluded. ChronoPack variants (rebuilt from the GRN step's components,
TF prior = priority_act share from tf_ranking.csv, target score = 1 - lfdr_bum):
  full                  TF prior x target dynamism (tie order: Wald statistic), = the package ranking
  TF prior only         TF prior for every candidate target (all targets of a TF tied)
  target dynamism only  1 - lfdr_bum (Wald tie order) for every TF (a gene's edges tied across TFs)
  TF prior x DE         TF prior x (1 - lfdr_bum) of conventional DE p-values (04a_de_target_scores.py)
  shuffled TF prior     TF prior values permuted among the universe TFs (N_SHUFFLE permutations)

Sensitivity: RELA/NFKB1 removed from the rewarded set (their edges stay in the ranking); minimum
reference targets 3 / 5 / 10; target universe restricted to genes scored by every method.

Writes output/tf_balanced_global/{tables,figures}/.
"""
from __future__ import annotations

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
from plot_regulator_prioritization_v2 import K562_EVAL_TFS, THP1_EVAL_TFS

OUT = fig5_paths.OUTPUT_DIR / "tf_balanced_global"
TABLES, FIGS = OUT / "tables", OUT / "figures"
DE_DIR = OUT / "de_scores"
COEXP = ["GENIE3", "GRNBOOST2", "PPCOR", "SINCERITIES"]
CONTEXT_TFS = {"K562": list(K562_EVAL_TFS), "THP1": list(THP1_EVAL_TFS)}
BENCHMARKS = ["TRRUST", "K562_ChIP"]
BUDGETS = [100, 300, 1000, 3000, 10000]
K_GRID = np.unique(np.concatenate([np.round(np.logspace(1, 5, 161)).astype(int), BUDGETS]))
LOGAUC_RANGE = (10, 10_000)
MIN_TARGETS = 3
MIN_TARGETS_SENS = [3, 5, 10]
EXCLUDED_DOMINANT = ["RELA", "NFKB1"]
MEANINGFUL_RECALL = 0.05  # a TF "recovers targets" at budget K if >= 5% of its known targets are in the top K
N_SHUFFLE = 200
N_BOOT = 2000
SEED = 20261010

VARIANT_LABELS = {
    "full": "ChronoPack (TF prior x dynamism)",
    "tf_only": "TF prior only",
    "target_only": "Target dynamism only",
    "tf_x_de": "TF prior x DE",
    "shuffled": "Shuffled TF prior (mean)",
}
COLORS = {"ChronoPack_K562": b.METHOD_COLORS.get("ChronoPack_K562", "#2f6f4e"), "ChronoPack_THP1": b.METHOD_COLORS.get("ChronoPack_THP1", "#1f5f8f"),
          **{m: b.METHOD_COLORS.get(m, "#888888") for m in COEXP}, "Random": "#9a9a9a"}
ABLATION_COLORS = {"full": "#111111", "tf_only": "#d95f02", "target_only": "#1b9e77", "tf_x_de": "#7570b3", "shuffled": "#e7298a"}


# ---------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------

def chronopack_components(context: str, valid: set[str]) -> tuple[pd.Series, pd.DataFrame, pd.DataFrame]:
    """TF prior share, tradeSeq target scores and DE target scores for one context."""
    tfs = pd.read_csv(fig5_paths.CHRONOPACK_OUTPUT_DIR / context / "tf_ranking.csv")  # Fig5_GRN_v4: step 01 ranking
    tfs["priority_act"] = pd.to_numeric(tfs["priority_act"], errors="coerce")
    tfs = tfs[np.isfinite(tfs["priority_act"]) & (tfs["priority_act"] > 0)].drop_duplicates("TF")
    prior = (tfs.set_index("TF")["priority_act"] / tfs["priority_act"].sum()).astype(float)

    def targets(path: Path) -> pd.DataFrame:
        t = pd.read_csv(path)
        t = t.assign(gene=t["gene"].astype(str), score=1.0 - pd.to_numeric(t["lfdr_bum"], errors="coerce"))
        t = t[np.isfinite(t["score"]) & (t["score"] > 0) & t["gene"].isin(valid)].drop_duplicates("gene")
        return t[["gene", "score", "waldStat_1"]].reset_index(drop=True)

    tradeseq = targets(fig5_paths.CHRONOPACK_OUTPUT_DIR / context / "target_scores.scores.csv")
    de = targets(DE_DIR / context / "target_scores.scores.csv")
    return prior, tradeseq, de


def outer_edges(prior: pd.Series, targets: pd.DataFrame, tfs: list[str], *, use_prior: bool, use_target: bool) -> pd.DataFrame:
    """Edges TF x target with key1 = (prior) x (target score); key2 = Wald/DE statistic when the target score is used."""
    p = prior.reindex(tfs).to_numpy(float)
    g = targets["score"].to_numpy(float)
    k1 = np.outer(p if use_prior else np.ones_like(p), g if use_target else np.ones_like(g)).ravel()
    k2 = np.tile(targets["waldStat_1"].to_numpy(float), len(tfs)) if use_target else np.zeros(len(tfs) * len(g))
    df = pd.DataFrame({"TF": np.repeat(tfs, len(g)), "target": np.tile(targets["gene"].to_numpy(), len(tfs)), "key1": k1, "key2": k2})
    return df[df["TF"] != df["target"]].reset_index(drop=True)


# ---------------------------------------------------------------------------------------------
# Recall engine (ties and unranked edges in expectation)
# ---------------------------------------------------------------------------------------------

class Ranking:
    """A method's ranking of candidate edges, reduced to tie groups."""

    def __init__(self, edges: pd.DataFrame, n_universe: int):
        cols = ["TF", "target", "key1", "key2"] + (["pos"] if "pos" in edges.columns else [])
        e = edges[cols].copy()
        e["key2"] = e["key2"].fillna(-np.inf)
        order = np.lexsort((-e["key2"].to_numpy(), -e["key1"].to_numpy()))
        e = e.iloc[order].reset_index(drop=True)
        k1, k2 = e["key1"].to_numpy(), e["key2"].to_numpy()
        new_group = np.ones(len(e), dtype=bool)
        new_group[1:] = (k1[1:] != k1[:-1]) | (k2[1:] != k2[:-1])
        gid = np.cumsum(new_group) - 1
        starts = np.flatnonzero(new_group)
        sizes = np.diff(np.append(starts, len(e)))
        self.edges = e.assign(group=gid)
        self.g_start, self.g_size = starts.astype(float), sizes.astype(float)
        self.n_emit, self.n_universe = len(e), n_universe

    def expected_counts(self, mask: np.ndarray, tf_of_row: np.ndarray, tfs: list[str], total_by_tf: dict[str, int], ks: np.ndarray) -> pd.DataFrame:
        """Expected number of masked edges per TF among the top K (rows: TF, columns: K).

        total_by_tf gives each TF's masked edges in the whole universe; those not emitted sit in the tail.
        """
        out = pd.DataFrame(0.0, index=tfs, columns=ks)
        sel = pd.DataFrame({"TF": tf_of_row[mask], "group": self.edges["group"].to_numpy()[mask]})
        counts = sel.groupby(["TF", "group"]).size()
        tail_start, tail_size = float(self.n_emit), float(self.n_universe - self.n_emit)
        for tf in tfs:
            frac_list, cnt_list = [], []
            if tf in counts.index.get_level_values(0):
                c = counts.loc[tf]
                s, n = self.g_start[c.index.to_numpy()], self.g_size[c.index.to_numpy()]
                frac = np.clip((ks[None, :] - s[:, None]) / n[:, None], 0, 1)
                out.loc[tf] = (c.to_numpy()[:, None] * frac).sum(axis=0)
                emitted = int(c.sum())
            else:
                emitted = 0
            rest = total_by_tf[tf] - emitted
            if rest > 0 and tail_size > 0:
                out.loc[tf] += rest * np.clip((ks - tail_start) / tail_size, 0, 1)
        return out

    def recall(self, positives: set[tuple[str, str]], pos_by_tf: dict[str, int], tfs: list[str], ks: np.ndarray) -> pd.DataFrame:
        if "pos" in self.edges.columns:  # precomputed membership (shuffle loop)
            mask = self.edges["pos"].to_numpy(bool)
        else:
            mask = pd.MultiIndex.from_frame(self.edges[["TF", "target"]]).isin(list(positives))
        tp = self.expected_counts(mask, self.edges["TF"].to_numpy(), tfs, pos_by_tf, ks)
        return tp.div(pd.Series(pos_by_tf).reindex(tfs), axis=0)

    def edge_share(self, tfs: list[str], n_targets: int, ks: np.ndarray) -> pd.DataFrame:
        """Expected number of each TF's edges among the top K (all candidate edges, annotated or not)."""
        mask = self.edges["TF"].isin(tfs).to_numpy()
        return self.expected_counts(mask, self.edges["TF"].to_numpy(), tfs, {t: n_targets for t in tfs}, ks)


def log_auc(curve: pd.Series) -> float:
    ks = curve.index.to_numpy(dtype=float)
    keep = (ks >= LOGAUC_RANGE[0]) & (ks <= LOGAUC_RANGE[1])
    x, y = np.log10(ks[keep]), curve.to_numpy()[keep]
    trapezoid = getattr(np, "trapezoid", None) or np.trapz
    return float(trapezoid(y, x) / (x[-1] - x[0]))


# ---------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------

def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    valid = b.load_valid_targets()
    gts = pub_v2.load_reference_ground_truths(valid)
    gt_tfs = set().union(*[set(g["TF"]) for g in gts.values()])

    coexp = {}
    for m in COEXP:
        e = compare.c.load_beeline_method(m, valid)
        coexp[m] = e[e["TF"].isin(gt_tfs)].rename(columns={"EdgeWeight": "key1"}).assign(key2=0.0)[["TF", "target", "key1", "key2"]]
    comps = {ctx: chronopack_components(ctx, valid) for ctx in CONTEXT_TFS}

    # Check: the rebuilt "full" ranking reproduces the package's GRN ranking (grn.edges.csv rank order).
    for ctx, (prior, ts, _de) in comps.items():
        pkg = pd.read_csv(fig5_paths.CHRONOPACK_OUTPUT_DIR / ctx / "grn.edges.csv", usecols=["regulator", "target", "rank"])
        tfs_chk = list(prior.index[:20])
        mine = Ranking(outer_edges(prior, ts, tfs_chk, use_prior=True, use_target=True), 1).edges
        ref = pkg[pkg["regulator"].isin(tfs_chk) & pkg["target"].isin(set(ts["gene"]))].sort_values("rank")
        same = (mine["TF"].to_numpy()[:20000] == ref["regulator"].to_numpy()[:20000]) & (mine["target"].to_numpy()[:20000] == ref["target"].to_numpy()[:20000])
        print(f"[check] {ctx}: rebuilt full ranking matches package order for {same.mean():.4%} of the first 20,000 edges (top-20 TFs)", flush=True)

    summary_rows, pertf_rows, share_rows, shuffle_rows, curve_rows = [], [], [], [], []
    for bench in BENCHMARKS:
        gt = gts[bench]
        avail = [set(e["TF"]) for e in coexp.values()] + [set(c[0].index) for c in comps.values()]
        universe_tfs = sorted(set(gt["TF"]).intersection(*avail))
        print(f"\n== {bench}: universe {len(universe_tfs)} TFs x {len(valid)} targets", flush=True)

        target_sets = {"all": sorted(valid)}
        shared = set(valid)
        for e in coexp.values():
            shared &= set(e["target"])
        for prior, ts, de in comps.values():
            shared &= set(ts["gene"]) & set(de["gene"])
        target_sets["shared"] = sorted(shared)
        print(f"   shared-target universe: {len(shared)} genes scored by every method", flush=True)

        for tset_name, targets_u in target_sets.items():
            tset = set(targets_u)
            n_universe = len(universe_tfs) * len(targets_u)
            gt_u = gt[gt["TF"].isin(universe_tfs) & gt["target"].isin(tset) & (gt["TF"] != gt["target"])].drop_duplicates(["TF", "target"])
            positives = set(zip(gt_u["TF"], gt_u["target"]))
            pos_by_tf_all = gt_u["TF"].value_counts().to_dict()

            def restrict(e: pd.DataFrame) -> pd.DataFrame:
                return e[e["TF"].isin(universe_tfs) & e["target"].isin(tset) & (e["TF"] != e["target"])]

            for ctx, (prior, ts, de) in comps.items():
                signaling = [t for t in CONTEXT_TFS[ctx] if t in universe_tfs]
                rewarded_sets = {}
                for mt in MIN_TARGETS_SENS:
                    tfs_ok = [t for t in signaling if pos_by_tf_all.get(t, 0) >= mt]
                    rewarded_sets[f"min{mt}"] = tfs_ok
                    rewarded_sets[f"min{mt}_noRELA_NFKB1"] = [t for t in tfs_ok if t not in EXCLUDED_DOMINANT]
                if tset_name == "shared":
                    rewarded_sets = {k: v for k, v in rewarded_sets.items() if k in (f"min{MIN_TARGETS}", f"min{MIN_TARGETS}_noRELA_NFKB1")}
                all_rewarded = sorted(set().union(*rewarded_sets.values()))
                print(f"   [{tset_name}] {ctx}: signaling TFs in universe {signaling}; targets "
                      + ", ".join(f"{t}={pos_by_tf_all.get(t, 0)}" for t in signaling), flush=True)

                cp_name = f"ChronoPack_{ctx}"
                rankings = {
                    f"{cp_name}|full": outer_edges(prior, ts, universe_tfs, use_prior=True, use_target=True),
                    f"{cp_name}|tf_only": outer_edges(prior, ts, universe_tfs, use_prior=True, use_target=False),
                    f"{cp_name}|target_only": outer_edges(prior, ts, universe_tfs, use_prior=False, use_target=True),
                    f"{cp_name}|tf_x_de": outer_edges(prior, de, universe_tfs, use_prior=True, use_target=True),
                    **{m: e for m, e in coexp.items()},
                }
                ks = K_GRID.astype(float)
                per_method_recall = {}
                for name, e in rankings.items():
                    r = Ranking(restrict(e), n_universe)
                    rec = r.recall(positives, {t: pos_by_tf_all.get(t, 0) for t in all_rewarded}, all_rewarded, ks)
                    per_method_recall[name] = rec
                    if tset_name == "all":
                        share = r.edge_share(universe_tfs, len(targets_u), np.array(BUDGETS, dtype=float))
                        for tf in universe_tfs:
                            share_rows.append({"benchmark": bench, "context": ctx, "method": name, "TF": tf,
                                               "signaling": tf in signaling, **{f"edges_top{k}": share.loc[tf, float(k)] for k in BUDGETS}})
                rand = pd.DataFrame(np.tile(ks / n_universe, (len(all_rewarded), 1)), index=all_rewarded, columns=ks)
                per_method_recall["Random"] = rand

                # Shuffled TF prior: permute prior values among universe TFs.
                shuf_curves = {k: [] for k in rewarded_sets}
                base_prior = prior.reindex(universe_tfs).to_numpy(float)
                base = restrict(outer_edges(prior, ts, universe_tfs, use_prior=True, use_target=True)).reset_index(drop=True)
                base = base.assign(pos=pd.MultiIndex.from_frame(base[["TF", "target"]]).isin(list(positives)))
                tf_code = pd.Categorical(base["TF"], categories=universe_tfs).codes
                gene_score = base["key1"].to_numpy() / base_prior[tf_code]  # target score of each edge
                for i in range(N_SHUFFLE if tset_name == "all" else 0):
                    perm = rng.permutation(base_prior)
                    e = base.assign(key1=perm[tf_code] * gene_score)
                    rec = Ranking(e, n_universe).recall(positives, {t: pos_by_tf_all.get(t, 0) for t in all_rewarded}, all_rewarded, ks)
                    for k, tfl in rewarded_sets.items():
                        if tfl:
                            shuf_curves[k].append(rec.loc[tfl].mean(axis=0))

                for set_name, tfl in rewarded_sets.items():
                    if not tfl:
                        continue
                    for name, rec in per_method_recall.items():
                        curve = rec.loc[tfl].mean(axis=0)
                        row = summarise(bench, ctx, tset_name, set_name, name, tfl, rec, curve, n_universe, rng)
                        summary_rows.append(row)
                        if set_name in (f"min{MIN_TARGETS}", f"min{MIN_TARGETS}_noRELA_NFKB1"):
                            for kk, v in curve.items():
                                curve_rows.append({"benchmark": bench, "target_universe": tset_name, "context": ctx, "rewarded": set_name, "method": name, "K": kk, "R_bal": v})
                    if not shuf_curves[set_name]:
                        continue
                    sc = np.vstack([c.to_numpy() for c in shuf_curves[set_name]])
                    mean_curve = pd.Series(sc.mean(axis=0), index=ks)
                    lo, hi = np.percentile(sc, [5, 95], axis=0)
                    aucs = np.array([log_auc(pd.Series(c, index=ks)) for c in sc])
                    full_auc = log_auc(per_method_recall[f"{cp_name}|full"].loc[tfl].mean(axis=0))
                    shuffle_rows.append({"benchmark": bench, "target_universe": tset_name, "context": ctx, "rewarded": set_name,
                                         "n_TFs": len(tfl), "logAUC_full": full_auc, "logAUC_shuffled_mean": aucs.mean(),
                                         "logAUC_shuffled_p5": np.percentile(aucs, 5), "logAUC_shuffled_p95": np.percentile(aucs, 95),
                                         "frac_shuffles_ge_full": float((aucs >= full_auc).mean()),
                                         **{f"R@{k}_shuffled_mean": float(mean_curve.loc[float(k)]) for k in BUDGETS}})
                    if set_name in (f"min{MIN_TARGETS}", f"min{MIN_TARGETS}_noRELA_NFKB1"):
                        for kk, v, l, h in zip(ks, mean_curve.to_numpy(), lo, hi):
                            curve_rows.append({"benchmark": bench, "target_universe": tset_name, "context": ctx, "rewarded": set_name,
                                               "method": f"{cp_name}|shuffled", "K": kk, "R_bal": v, "lo": l, "hi": h})

                if tset_name == "all":
                    tfl = rewarded_sets[f"min{MIN_TARGETS}"]
                    for name, rec in per_method_recall.items():
                        for tf in tfl:
                            emitted = rankings[name] if name in rankings else None
                            n_emit_pos = np.nan
                            if emitted is not None:
                                e = restrict(emitted)
                                n_emit_pos = int(pd.MultiIndex.from_frame(e[e["TF"] == tf][["TF", "target"]]).isin(list(positives)).sum())
                            pertf_rows.append({"benchmark": bench, "context": ctx, "method": name, "TF": tf,
                                               "n_ref_targets": pos_by_tf_all.get(tf, 0), "n_ref_targets_ranked_by_method": n_emit_pos,
                                               **{f"recall@{k}": float(rec.loc[tf, float(k)]) for k in BUDGETS},
                                               **{f"contribution@{k}": float(rec.loc[tf, float(k)]) / len(tfl) for k in BUDGETS},
                                               "logAUC": log_auc(rec.loc[tf])})

    summary = pd.DataFrame(summary_rows)
    pertf = pd.DataFrame(pertf_rows)
    shares = pd.DataFrame(share_rows)
    shuffles = pd.DataFrame(shuffle_rows)
    curves = pd.DataFrame(curve_rows)
    summary.to_csv(TABLES / "summary.csv", index=False)
    pertf.to_csv(TABLES / "per_tf_recall.csv", index=False)
    shares.to_csv(TABLES / "edge_share_by_tf.csv", index=False)
    shuffles.to_csv(TABLES / "shuffled_tf_prior.csv", index=False)
    curves.to_csv(TABLES / "curves.csv", index=False)
    make_figures(summary, pertf, curves)
    print(f"\nWrote {TABLES} and {FIGS}")


def summarise(bench, ctx, tset_name, set_name, name, tfl, rec, curve, n_universe, rng) -> dict:
    sub = rec.loc[tfl]
    tf_auc = np.array([log_auc(sub.loc[t]) for t in tfl])
    boot_idx = rng.integers(0, len(tfl), size=(N_BOOT, len(tfl)))
    boot_auc = tf_auc[boot_idx].mean(axis=1)
    row = {"benchmark": bench, "context": ctx, "target_universe": tset_name, "rewarded": set_name, "method": name,
           "n_TFs": len(tfl), "TFs": ",".join(tfl), "n_universe_edges": n_universe,
           "logAUC": log_auc(curve), "logAUC_ci_lo": np.percentile(boot_auc, 2.5), "logAUC_ci_hi": np.percentile(boot_auc, 97.5)}
    rand_auc = log_auc(pd.Series(curve.index.to_numpy(float) / n_universe, index=curve.index))
    row["logAUC_vs_random"] = row["logAUC"] / rand_auc
    for k in BUDGETS:
        vals = sub[float(k)].to_numpy()
        boot = vals[boot_idx].mean(axis=1)
        row[f"R@{k}"] = float(curve.loc[float(k)])
        row[f"R@{k}_ci_lo"], row[f"R@{k}_ci_hi"] = np.percentile(boot, [2.5, 97.5])
        row[f"n_TFs_recall>={MEANINGFUL_RECALL:g}@{k}"] = int((vals >= MEANINGFUL_RECALL).sum())
    return row


# ---------------------------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------------------------

def _label(name: str) -> str:
    if "|" in name:
        return VARIANT_LABELS[name.split("|")[1]]
    return name


def make_figures(summary: pd.DataFrame, pertf: pd.DataFrame, curves: pd.DataFrame) -> None:
    plt.rcParams.update({"pdf.fonttype": 42, "font.size": 9})
    main_set = f"min{MIN_TARGETS}"

    def grid(rewarded: str, tset: str, methods_fn, stem: str, title: str) -> None:
        fig, axes = plt.subplots(len(BENCHMARKS), 2, figsize=(10, 7.2), sharex=True)
        for i, bench in enumerate(BENCHMARKS):
            for j, ctx in enumerate(CONTEXT_TFS):
                ax = axes[i, j]
                c = curves[(curves.benchmark == bench) & (curves.context == ctx) & (curves.rewarded == rewarded) & (curves.target_universe == tset)]
                if c.empty:
                    ax.set_visible(False)
                    continue
                for name, color, lw, ls in methods_fn(ctx):
                    d = c[c.method == name].sort_values("K")
                    if d.empty:
                        continue
                    ax.plot(d.K, d.R_bal, color=color, lw=lw, ls=ls, label=_label(name) if name != f"ChronoPack_{ctx}|full" else f"ChronoPack ({ctx})")
                    if "lo" in d and d["lo"].notna().any():
                        ax.fill_between(d.K, d["lo"], d["hi"], color=color, alpha=0.15, lw=0)
                tfs = summary[(summary.benchmark == bench) & (summary.context == ctx) & (summary.rewarded == rewarded) & (summary.target_universe == tset)]["TFs"].iloc[0]
                ax.set_xscale("log")
                ax.set_title(f"{bench} | {ctx} context\nrewarded: {tfs}", fontsize=8.5, loc="left")
                ax.grid(True, color="#e8e8e8", lw=0.6)
                ax.spines[["top", "right"]].set_visible(False)
                if i == len(BENCHMARKS) - 1:
                    ax.set_xlabel("Global edge budget K (top-ranked edges)")
                if j == 0:
                    ax.set_ylabel("TF-balanced recall")
                ax.legend(fontsize=7, frameon=False, loc="upper left")
        fig.suptitle(title, fontsize=11)
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(FIGS / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
        plt.close(fig)

    def methods_main(ctx):
        return [(f"ChronoPack_{ctx}|full", COLORS[f"ChronoPack_{ctx}"], 2.2, "-")] + [(m, COLORS[m], 1.3, "-") for m in COEXP] + [("Random", COLORS["Random"], 1.2, "--")]

    def methods_ablation(ctx):
        return [(f"ChronoPack_{ctx}|{v}", ABLATION_COLORS[v], 2.0 if v == "full" else 1.4, "-") for v in VARIANT_LABELS] + [("Random", COLORS["Random"], 1.2, "--")]

    grid(main_set, "all", methods_main, "1_tf_balanced_recall_methods", "TF-balanced global recall of signaling-TF targets (all 45/28 universe TFs ranked)")
    grid(f"{main_set}_noRELA_NFKB1", "all", methods_main, "4a_tf_balanced_recall_methods_noRELA_NFKB1", "Excluding RELA/NFKB1 from the rewarded set (their edges stay ranked)")
    grid(main_set, "shared", methods_main, "4b_tf_balanced_recall_methods_shared_targets", "Target universe restricted to genes scored by every method")
    grid(main_set, "all", methods_ablation, "5_chronopack_ablation", "ChronoPack score ablations (shuffled: mean and 5-95% band)")
    grid(f"{main_set}_noRELA_NFKB1", "all", methods_ablation, "5b_chronopack_ablation_noRELA_NFKB1", "ChronoPack ablations, RELA/NFKB1 not rewarded")

    # Per-TF recall heatmaps at two budgets.
    for bench in BENCHMARKS:
        fig, axes = plt.subplots(1, 2 * len(CONTEXT_TFS), figsize=(16, 4.2))
        col = 0
        for ctx in CONTEXT_TFS:
            p = pertf[(pertf.benchmark == bench) & (pertf.context == ctx)]
            if p.empty:
                axes[col].set_visible(False); axes[col + 1].set_visible(False); col += 2
                continue
            for k in (1000, 10000):
                ax = axes[col]
                mat = p.pivot(index="method", columns="TF", values=f"recall@{k}")
                order = [m for m in [f"ChronoPack_{ctx}|{v}" for v in ("full", "tf_only", "target_only", "tf_x_de")] + COEXP + ["Random"] if m in mat.index]
                mat = mat.loc[order]
                ntg = p.drop_duplicates("TF").set_index("TF")["n_ref_targets"]
                im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="viridis", vmin=0, vmax=max(0.05, float(np.nanmax(mat.to_numpy()))))
                ax.set_xticks(range(mat.shape[1]))
                ax.set_xticklabels([f"{t}\n(n={ntg[t]})" for t in mat.columns], fontsize=7)
                ax.set_yticks(range(mat.shape[0]))
                ax.set_yticklabels([_label(m) for m in mat.index] if col == 0 else [], fontsize=7)
                for (r, c_), v in np.ndenumerate(mat.to_numpy()):
                    ax.text(c_, r, f"{v:.2f}", ha="center", va="center", fontsize=6, color="white" if v < 0.6 * im.get_clim()[1] else "black")
                ax.set_title(f"{ctx} context, recall@{k}", fontsize=9)
                col += 1
        fig.suptitle(f"{bench}: per-TF recall of known targets among the top K global edges", fontsize=10)
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(FIGS / f"3_per_tf_recall_{bench}.{ext}", dpi=200, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
