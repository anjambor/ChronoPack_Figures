#!/usr/bin/env python
"""Perturb-seq as a functional reference: score the with_prior GRNs against Replogle (K562) and Yao (THP1) TF
perturbation responses (user request 2026-11-11), instead of using the perturbation rankings as comparator methods.

References (perturbed TF vs non-targeting control, Seurat FindMarkers over all measured genes; the TF's own gene is
excluded):
  Replogle2022 K562  input/capsule/external/.../deg_tables (+ FOS from input/capsule/perturb_de/replogle_fos): 11 TFs
  Yao2023 THP1 LPS   causalbench Yao2023 TFs_vs_non_targeting_DE_all.csv (+ IRF3 from input/capsule/perturb_de/yao_irf3): 12 TFs
Methods: the with_prior networks without Replogle / Yao (ChronoPack K562 / THP1, BEELINE, pySCENIC, CellOracle,
Inferelator), from the step 02 ranked-edge tables (the edges 5b scores). A TF a method has no predictions for is left
out for that method, as in 5b.
Gene universes (user choice: both):
  shared  the 2,077-gene shared universe of the *_v2 figures (tables/shared_universe_v2), within the genes measured
          in the reference
  full    Fig5's 19,338 protein-coding genes within the genes measured in the reference
Metrics (user choice: rank-based primary, FDR secondary):
  rank-based, no threshold on the response: the response score of a gene is -log10 P of its DE test (direction-
  agnostic: most methods do not predict a sign). For each method and TF, the predicted target set is the method's
  top K targets within the universe (K = 100; 50 and 200 in the table; all predicted targets if fewer, at least 10).
    response AUROC  P(response of a predicted target > response of another universe gene) (Mann-Whitney; 0.5 = random)
    GSEA            weighted (p = 1) GSEA enrichment score (ES) of the predicted set in the genes ranked by response,
                    against N_PERM random sets of the same size from the universe: ES z-score (plotted; 0 = random),
                    NES (ES / mean null ES of the same sign; random sets give |NES| ~ 1, so NES is tabled only) and the
                    permutation P (same-sign nulls), BH within each panel
  BH FDR < 0.05: positives = genes with Benjamini-Hochberg q < 0.05 (per TF over all genes measured in the reference;
  the tables' own adjusted P, Seurat p_val_adj, is Bonferroni and was used before 2026-11-11); per-TF AUPRCr and EPR
  with Fig5's evaluate_one (random tail for unranked genes), only for TFs with >= MIN_POSITIVES positives.
  Threshold-free (user request): soft labels = 1 - local FDR (beta-uniform mixture fit to each TF's P values, Pounds &
  Morris 2003; lfdr = pi0 / f(p)), giving soft-label ROC curves (TPR / FPR = cumulative responder / non-responder mass);
  and a threshold sweep (AUROC with the top N responders by P as positives, N = SWEEP_N).
Curves (user request): GSEA running-enrichment curves of each method's top-K set per TF over a random-set band, and top-k response curves
(mean response percentile of a method's top k targets, k = 1..500, averaged over TFs; 0.5 = random).
scMINER Fig. 4-style panels (user request): ROC curves with BH FDR < 0.05 DEGs as labels (global = all TF-gene pairs
of a method's TFs pooled on edge weight; per TF for TFs with >= MIN_POSITIVES DEGs), a summary of per-method medians
over TFs (TF-bootstrap 95% CI, TF dots, paired Wilcoxon vs the in-context ChronoPack) and a hit rug under the GSEA
curves. Signed (positive / negative target) GSEA is not possible: ChronoPack's signs are per TF, not per edge.
``--curves-only`` redraws curves, ROC and summary panels from the existing metric tables.
These responses include indirect effects of the TF, so they measure functional influence, not direct binding.
ChronoPack ranks targets in the same order for every TF (edge weight = TF priority share x target score), so its
per-TF scores test its target score against each TF's response.

Outputs: output/with_prior/benchmark_outputs/figures/perturbseq_reference_v2/ and tables/perturbseq_reference_v2/.
Run after the with_prior steps 02 (full) and 02d (shared): FIG5_METHOD_SET=with_prior python scripts/02e_perturbseq_reference.py
"""
from __future__ import annotations

import glob
import json
import os
import pickle
import sys
import zlib
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
for p in (SCRIPTS, SCRIPTS / "capsule_lib" / "01_generate_predictions_and_metrics"):
    sys.path.insert(0, str(p))

if os.environ.get("FIG5_METHOD_SET") != "with_prior" or os.environ.get("FIG5_UNIVERSE", "full") != "full":
    sys.exit("run with FIG5_METHOD_SET=with_prior (and FIG5_UNIVERSE unset)")

import fig5_paths  # noqa: E402

os.environ.setdefault("MPLCONFIGDIR", str(fig5_paths.OUTPUT_DIR / ".mplconfig"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from scipy.stats import mannwhitneyu, rankdata  # noqa: E402

import benchmark_grn_systems as b  # noqa: E402

FIG_OUT = fig5_paths.FIGURE_DIR / "perturbseq_reference_v2"
TAB_OUT = fig5_paths.TABLE_DIR / "perturbseq_reference_v2"
EDGES = {
    "full": fig5_paths.TABLE_DIR / "regulator_prioritization_method_specific_ranked_edges.csv",
    "shared": fig5_paths.TABLE_DIR / "shared_universe_v2" / "regulator_prioritization_method_specific_ranked_edges.csv",
}
SHARED_GENES = fig5_paths.TABLE_DIR / "shared_universe_v2" / "shared_universe_genes.txt"
REPLOGLE_FOS = fig5_paths.CAPSULE_INPUT_DIR / "perturb_de" / "replogle_fos" / "deg_tables" / "ENSG00000170345.findmarkers.all.txt"
YAO_IRF3 = fig5_paths.CAPSULE_INPUT_DIR / "perturb_de" / "yao_irf3" / "TFs_vs_non_targeting_DE_IRF3.csv"

METHODS = [m for m in b.METHOD_ORDER if m not in fig5_paths.PERTURB_METHODS]
THP1_METHODS = {"ChronoPack_THP1"}
REFERENCES = {"Replogle_K562": "Replogle2022 | K562 Perturb-seq", "Yao_THP1": "Yao2023 | THP1 LPS Perturb-seq"}
REFERENCE_CONTEXT = {"Replogle_K562": "K562", "Yao_THP1": "THP1"}
UNIVERSE_LABELS = {"shared": "shared 2,077-gene universe", "full": "all protein-coding genes"}
K_MAIN, K_ALL, MIN_SET = 100, (50, 100, 200), 10
N_PERM, SEED = 1000, 777
FDR, MIN_POSITIVES = 0.05, 5
SHORT = {"ChronoPack_K562": "Chrono", "ChronoPack_THP1": "Chrono", **fig5_paths.PRIOR_LABELS}
CMAP = LinearSegmentedColormap.from_list("div", ["#2b6a99", "#f2f2f2", "#c0504d"])


def load_references() -> dict[str, pd.DataFrame]:
    """reference -> (TF, target, p_val, log2FC, p_bonf, q_bh, lfdr, soft, pi0); one row per measured gene, the TF's own
    gene removed. p_bonf is the tables' Seurat p_val_adj (Bonferroni); q_bh, lfdr and pi0 are per TF over its genes."""
    with open(b.EN2SYM_PATH, "rb") as f:
        en2sym = pickle.load(f)
    rows = []
    files = [(Path(p), None) for p in sorted(glob.glob(str(b.REPLOGLIE_DEG_DIR / "*")))] + [(REPLOGLE_FOS, "FOS")]
    for path, tf in files:
        tf = tf or en2sym.get(path.name.split(".")[0])
        t = pd.read_csv(path, sep="\t")
        rows.append(pd.DataFrame({"TF": tf, "target": t.iloc[:, 0].astype(str).map(en2sym), "p_val": t["p_val"],
                                  "log2FC": t["avg_log2FC"], "p_bonf": t["p_val_adj"]}))
    rep = pd.concat(rows, ignore_index=True)
    yao = pd.concat([pd.read_csv(b.YAO2023_DEGS_PATH), pd.read_csv(YAO_IRF3)], ignore_index=True)
    yao = pd.DataFrame({"TF": yao["TF"], "target": yao["Gene name"], "p_val": yao["P-value"],
                        "log2FC": yao["Avg_log2FC"], "p_bonf": yao["P_val_adj"]})
    out = {}
    for name, df in {"Replogle_K562": rep, "Yao_THP1": yao}.items():
        df = df.dropna(subset=["TF", "target", "p_val"])
        df = df[df["TF"] != df["target"]]
        df = df.sort_values("p_val", kind="mergesort").drop_duplicates(["TF", "target"]).reset_index(drop=True)
        df["score"] = -np.log10(df["p_val"].clip(lower=np.nextafter(0, 1)))
        df["q_bh"] = df.groupby("TF")["p_val"].transform(lambda v: bh(v).to_numpy())
        lf = {tf: bum_lfdr(g["p_val"].to_numpy()) for tf, g in df.groupby("TF")}
        df["lfdr"] = np.nan
        for tf, (lfdr, _) in lf.items():
            df.loc[df["TF"] == tf, "lfdr"] = lfdr
        df["soft"] = 1 - df["lfdr"]
        df["pi0"] = df["TF"].map({tf: v[1] for tf, v in lf.items()})
        out[name] = df
    return out


def load_edges(universe: str) -> dict[tuple[str, str], list[str]]:
    """(method, TF) -> targets in the method's rank order (step 02 ranked-edge table)."""
    e = pd.read_csv(EDGES[universe], usecols=["method", "TF", "target", "EdgeWeight"])
    e = e[e["method"].isin(METHODS)]
    e = e.sort_values("EdgeWeight", ascending=False, kind="mergesort")
    return {key: g["target"].tolist() for key, g in e.groupby(["method", "TF"], sort=False)}


def gsea_es(order_weights: np.ndarray, hits: np.ndarray) -> np.ndarray:
    """Weighted (p = 1) GSEA enrichment score; hits is (n_sets, n_genes) boolean in response order."""
    hit_w = hits * order_weights
    hit_sum = hit_w.sum(axis=1, keepdims=True)
    miss = (~hits).astype(float)
    run = np.cumsum(hit_w / np.where(hit_sum > 0, hit_sum, 1) - miss / miss.sum(axis=1, keepdims=True), axis=1)
    mx, mn = run.max(axis=1), run.min(axis=1)
    return np.where(np.abs(mx) >= np.abs(mn), mx, mn)


class Null:
    """Cached GSEA nulls: random sets of a given size from one (universe, reference, TF) ranking."""

    def __init__(self) -> None:
        self.cache: dict[tuple, np.ndarray] = {}

    def get(self, key: tuple, weights: np.ndarray, size: int) -> np.ndarray:
        if (key, size) not in self.cache:
            rng = np.random.default_rng([SEED, size, zlib.crc32(repr(key).encode())])
            idx = np.argpartition(rng.random((N_PERM, len(weights))), size - 1, axis=1)[:, :size]
            hits = np.zeros((N_PERM, len(weights)), dtype=bool)
            np.put_along_axis(hits, idx, True, axis=1)
            self.cache[(key, size)] = gsea_es(weights, hits)
        return self.cache[(key, size)]


def bum_lfdr(p: np.ndarray) -> tuple[np.ndarray, float]:
    """Local FDR from P values with a beta-uniform mixture (BUM; Pounds & Morris 2003), fitted by maximum likelihood:
    f(p) = lam + (1 - lam) a p^(a - 1), 0 < a < 1; pi0 = lam + (1 - lam) a (the BUM upper bound); lfdr = min(1, pi0 / f)."""
    from scipy.optimize import minimize

    logp = np.log(np.clip(np.asarray(p, float), 1e-300, 1.0))

    def logf(theta):
        lam, a = theta
        return np.logaddexp(np.log(lam), np.log1p(-lam) + np.log(a) + (a - 1) * logp)

    eps = 1e-6
    fit = minimize(lambda t: -logf(t).sum(), x0=[0.7, 0.5], method="L-BFGS-B", bounds=[(eps, 1 - eps), (eps, 1 - eps)])
    lam, a = fit.x
    pi0 = float(min(1.0, lam + (1 - lam) * a))
    return np.minimum(1.0, pi0 / np.exp(logf(fit.x))), pi0


def bh(p: pd.Series) -> pd.Series:
    p = p.astype(float)
    ok = p.notna()
    q = pd.Series(np.nan, index=p.index)
    if ok.any():
        v = p[ok].to_numpy()
        order = np.argsort(v)
        adj = np.minimum.accumulate((v[order] * len(v) / np.arange(1, len(v) + 1))[::-1])[::-1]
        q.loc[p[ok].index[order]] = np.minimum(adj, 1)
    return q


def rank_metrics(refs, edges, universes) -> pd.DataFrame:
    null = Null()
    rows = []
    for uname, valid in universes.items():
        for rname, ref in refs.items():
            for tf, r in ref.groupby("TF", sort=False):
                r = r[r["target"].isin(valid)]
                # response order (ties broken by gene name, deterministic)
                r = r.sort_values(["score", "target"], ascending=[False, True], kind="mergesort").reset_index(drop=True)
                genes = r["target"].to_numpy()
                weights = r["score"].to_numpy()
                pos = {g: i for i, g in enumerate(genes)}
                ranks = rankdata(weights)  # for the Mann-Whitney AUROC (ties averaged)
                for method in METHODS:
                    pred = [g for g in edges[uname].get((method, tf), []) if g in pos]
                    for k in K_ALL:
                        top = pred[:k]
                        if len(top) < MIN_SET:
                            continue
                        idx = np.array([pos[g] for g in top])
                        n1, n0 = len(idx), len(genes) - len(idx)
                        auroc = (ranks[idx].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
                        hits = np.zeros((1, len(genes)), dtype=bool)
                        hits[0, idx] = True
                        es = float(gsea_es(weights, hits)[0])
                        nul = null.get((uname, rname, tf), weights, len(idx))
                        same = nul[nul >= 0] if es >= 0 else nul[nul < 0]
                        nes = es / abs(same.mean()) if len(same) else np.nan
                        z = (es - nul.mean()) / nul.std()
                        p = (1 + (same >= es).sum() if es >= 0 else 1 + (same <= es).sum()) / (1 + len(same))
                        rows.append(dict(universe=uname, reference=rname, TF=tf, method=method, K=k, set_size=n1,
                                         n_universe=len(genes), n_predicted_in_universe=len(pred),
                                         response_AUROC=auroc, ES=es, ES_z=z, NES=nes, perm_p=p,
                                         context_match=(method in THP1_METHODS) == (REFERENCE_CONTEXT[rname] == "THP1")))
    m = pd.DataFrame(rows)
    m["perm_q"] = m.groupby(["universe", "reference", "K"])["perm_p"].transform(bh)
    return m


def fdr_metrics(refs, edges_df, universes) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, counts = [], []
    for uname, valid in universes.items():
        for rname, ref in refs.items():
            for tf, r in ref.groupby("TF", sort=False):
                u = set(r["target"]) & valid
                gt = r[(r["q_bh"] < FDR) & r["target"].isin(u)][["TF", "target"]]
                counts.append(dict(universe=uname, reference=rname, TF=tf, n_universe=len(u), n_positives=len(gt),
                                   n_positives_all_measured=int((r["q_bh"] < FDR).sum()),
                                   n_positives_bonferroni=int(((r["p_bonf"] < FDR) & r["target"].isin(u)).sum()),
                                   soft_positives=float(r.loc[r["target"].isin(u), "soft"].sum()), pi0=float(r["pi0"].iloc[0])))
                if len(gt) < MIN_POSITIVES:
                    continue
                for method in METHODS:
                    e = edges_df[uname]
                    e = e[(e["method"] == method) & (e["TF"] == tf)]
                    if e.empty:
                        continue
                    row, _ = b.evaluate_one(e[["TF", "target", "EdgeWeight"]], gt, {tf}, u, curve_points=2)
                    row.update(universe=uname, reference=rname, TF=tf, method=method)
                    rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(counts)


def heatmap_panel(ax, mat: pd.DataFrame, norm, fmt, marks: pd.DataFrame | None, title: str) -> None:
    rows = list(mat.index)
    ax.imshow(np.ma.masked_invalid(mat.to_numpy(float)), cmap=CMAP, norm=norm, aspect="auto")
    for i, method in enumerate(rows):
        for j, tf in enumerate(mat.columns):
            v = mat.iat[i, j]
            if pd.isna(v):
                continue
            star = "*" if marks is not None and tf in marks.columns and marks.at[method, tf] else ""
            strong = abs(norm(v) - 0.5) > 0.32
            ax.text(j, i, fmt(v) + star, ha="center", va="center", fontsize=7.4,
                    color="white" if strong else "#222222", fontweight="bold" if tf == "median" else None)
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([SHORT.get(m, m) for m in rows], fontsize=9)
    ax.set_title(title, fontsize=10, pad=6)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    if "median" in mat.columns:
        ax.axvline(len(mat.columns) - 1.5, color="#333333", lw=1.0)
    k562 = [i for i, m in enumerate(rows) if m not in THP1_METHODS]
    if k562 and len(k562) < len(rows):
        ax.axhline(max(k562) + 0.5, color="#333333", lw=1.0)
    for i in range(len(mat.columns) + 1):
        ax.axvline(i - 0.5, color="white", lw=0.8)
    for i in range(len(rows) + 1):
        ax.axhline(i - 0.5, color="white", lw=0.8)


def context_brackets(ax, rows: list[str]) -> None:
    k562 = [i for i, m in enumerate(rows) if m not in THP1_METHODS]
    thp1 = [i for i, m in enumerate(rows) if m in THP1_METHODS]
    for idx, label in ((k562, "K562 networks"), (thp1, "THP1")):
        if idx:  # fixed offset in points left of the axis, clear of the method labels
            ax.annotate(label, xy=(0, float(np.mean(idx))), xycoords=("axes fraction", "data"), xytext=(-78, 0),
                        textcoords="offset points", rotation=90, ha="center", va="center", fontsize=8.5,
                        color="#555555", annotation_clip=False)


def plot_rank(m: pd.DataFrame, universe: str, metric: str) -> None:
    sub = m[(m["universe"] == universe) & (m["K"] == K_MAIN)]
    if metric == "response_AUROC":
        span = max(0.05, float(np.nanmax(np.abs(sub[metric] - 0.5))))
        norm, fmt, label = TwoSlopeNorm(0.5, 0.5 - span, 0.5 + span), (lambda v: f"{v:.2f}"), "response AUROC (0.5 = random)"
    else:
        span = max(1.0, float(np.nanmax(np.abs(sub[metric]))))
        norm, fmt, label = TwoSlopeNorm(0.0, -span, span), (lambda v: f"{v:.1f}"), "GSEA ES z-score vs random sets (0 = random)"
    widths = [sub[sub["reference"] == r]["TF"].nunique() + 1 for r in REFERENCES]
    rows = [x for x in METHODS if x in set(sub["method"])]  # one row list for both panels
    fig, axes = plt.subplots(1, 2, figsize=(2.2 + 0.5 * sum(widths), 4.2), gridspec_kw={"width_ratios": widths})
    for ax, ref in zip(axes, REFERENCES):
        s = sub[sub["reference"] == ref]
        tfs = [t for t in dict.fromkeys(s.sort_values("TF")["TF"])]
        mat = s.pivot(index="method", columns="TF", values=metric).reindex(index=rows, columns=tfs)
        mat["median"] = mat.median(axis=1)
        marks = s.pivot(index="method", columns="TF", values="perm_q").reindex(index=mat.index, columns=tfs) < 0.05
        heatmap_panel(ax, mat, norm, fmt, marks, REFERENCES[ref])
        if ax is axes[1]:
            ax.set_yticklabels([])
    context_brackets(axes[0], list(mat.index))
    cax = fig.add_axes([0.93, 0.25, 0.012, 0.5])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label(label, fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    fig.suptitle(f"Perturbation response of each TF's top {K_MAIN} predicted targets ({UNIVERSE_LABELS[universe]})",
                 fontsize=10.5, y=0.995)
    fig.text(0.5, 0.005, "* GSEA permutation FDR < 0.05 (BH within panel). Blank: no predictions for that TF. "
             "Responses include indirect effects.", ha="center", fontsize=7.5, color="#555555")
    fig.subplots_adjust(left=0.1, right=0.9, top=0.86, bottom=0.2, wspace=0.05)
    stem = FIG_OUT / f"perturbseq_rank_{metric}_{universe}"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_fdr(f: pd.DataFrame, counts: pd.DataFrame, universe: str) -> None:
    sub = f[f["universe"] == universe] if len(f) else f
    refs = [r for r in REFERENCES if len(sub) and (sub["reference"] == r).any()]
    if not refs:
        return
    vals = np.log2(sub["AUPRC_ratio"].clip(lower=1e-3))
    span = max(1.0, float(np.nanmax(np.abs(vals))))
    norm = TwoSlopeNorm(0.0, -span, span)
    widths = [max(2, sub[sub["reference"] == r]["TF"].nunique() + 1) for r in refs]
    rows = [x for x in METHODS if x in set(sub["method"])]  # one row list for all panels
    fig, axes = plt.subplots(1, len(refs), figsize=(2.6 + 0.55 * sum(widths), 4.2), gridspec_kw={"width_ratios": widths}, squeeze=False)
    for ax, ref in zip(axes[0], refs):
        s = sub[sub["reference"] == ref]
        c = counts[(counts["universe"] == universe) & (counts["reference"] == ref)].set_index("TF")["n_positives"]
        tfs = sorted(s["TF"].unique(), key=lambda t: -c.get(t, 0))
        mat = s.pivot(index="method", columns="TF", values="AUPRC_ratio").reindex(index=rows, columns=tfs)
        if len(tfs) > 1:
            mat["median"] = mat.median(axis=1)
        heatmap_panel(ax, np.log2(mat), norm, lambda v: f"{2 ** v:.2f}", None, REFERENCES[ref])
        ax.set_xticklabels([f"{t} ({int(c[t])})" if t in c else t for t in mat.columns], rotation=45, ha="right", fontsize=8.5)
        if ax is not axes[0][0]:
            ax.set_yticklabels([])
    context_brackets(axes[0][0], list(mat.index))
    cax = fig.add_axes([0.93, 0.25, 0.012, 0.5])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    cb = fig.colorbar(sm, cax=cax)
    ticks = [t for t in range(-6, 7) if -span <= t <= span]
    cb.set_ticks(ticks)
    cb.set_ticklabels([f"{2 ** t:g}" for t in ticks])
    cb.set_label("AUPRCr (1 = random)", fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    fig.suptitle(f"BH FDR < 0.05 perturbation DEGs as positives ({UNIVERSE_LABELS[universe]}; TFs with >= {MIN_POSITIVES} positives)",
                 fontsize=10.5, y=0.995)
    fig.text(0.5, 0.005, "Column labels: TF (number of positives in the universe). Blank: no predictions for that TF.",
             ha="center", fontsize=7.5, color="#555555")
    fig.subplots_adjust(left=0.12, right=0.9, top=0.86, bottom=0.22, wspace=0.05)
    stem = FIG_OUT / f"perturbseq_fdr05_AUPRCr_{universe}"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


CURVE_KMAX, CURVE_POINTS, N_BAND = 500, 400, 200
LEGEND = {"ChronoPack_K562": "ChronoPack | K562", "ChronoPack_THP1": "ChronoPack | THP1", **fig5_paths.PRIOR_LABELS}


def response_curves(refs, edges, universes) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per (universe, reference, TF, method):
    - GSEA running enrichment of the top-K_MAIN set along the genes ranked by response (CURVE_POINTS grid + extremes);
    - top-k response curve: mean response percentile (0 = weakest, 1 = strongest; ties averaged, 0.5 = random) of the
      method's top k targets, k = 1..min(CURVE_KMAX, predicted targets in the universe)."""
    run_rows, topk_rows, band_rows, hit_rows = [], [], [], []
    for uname, valid in universes.items():
        for rname, ref in refs.items():
            for tf, r in ref.groupby("TF", sort=False):
                r = r[r["target"].isin(valid)]
                r = r.sort_values(["score", "target"], ascending=[False, True], kind="mergesort").reset_index(drop=True)
                genes, weights = r["target"].to_numpy(), r["score"].to_numpy()
                n = len(genes)
                pos = {g: i for i, g in enumerate(genes)}
                pct = (rankdata(weights) - 1) / (n - 1)
                grid = np.unique(np.linspace(0, n - 1, CURVE_POINTS).astype(int))
                # random-set reference: running ES of N_BAND random sets of size K_MAIN (weighted GSEA on -log10 P
                # gives random sets clearly positive curves, so the band is the baseline to read the methods against)
                rng = np.random.default_rng([SEED, zlib.crc32(repr((uname, rname, tf)).encode())])
                hits = np.zeros((N_BAND, n), dtype=bool)
                np.put_along_axis(hits, np.argpartition(rng.random((N_BAND, n)), K_MAIN - 1, axis=1)[:, :K_MAIN], True, axis=1)
                hw = hits * weights
                nul = np.cumsum(hw / hw.sum(axis=1, keepdims=True) - (~hits) / (~hits).sum(axis=1, keepdims=True), axis=1)[:, grid]
                q05, q50, q95 = np.percentile(nul, [5, 50, 95], axis=0)
                band_rows.append(pd.DataFrame({"universe": uname, "reference": rname, "TF": tf, "gene_rank": grid + 1,
                                               "rank_fraction": (grid + 1) / n, "q05": q05, "q50": q50, "q95": q95}))
                for method in METHODS:
                    pred = [g for g in edges[uname].get((method, tf), []) if g in pos]
                    if len(pred) < MIN_SET:
                        continue
                    idx = np.array([pos[g] for g in pred[:K_MAIN]])
                    hit = np.zeros(n, dtype=bool)
                    hit[idx] = True
                    hw = hit * weights
                    run = np.cumsum(hw / hw.sum() - (~hit) / (~hit).sum())
                    keep = np.unique(np.r_[np.linspace(0, n - 1, CURVE_POINTS).astype(int), run.argmax(), run.argmin()])
                    run_rows.append(pd.DataFrame({"universe": uname, "reference": rname, "TF": tf, "method": method,
                                                  "gene_rank": keep + 1, "rank_fraction": (keep + 1) / n, "running_ES": run[keep]}))
                    hit_rows.append(pd.DataFrame({"universe": uname, "reference": rname, "TF": tf, "method": method,
                                                  "target": genes[np.sort(idx)], "gene_rank": np.sort(idx) + 1,
                                                  "rank_fraction": (np.sort(idx) + 1) / n}))
                    kk = np.arange(1, min(CURVE_KMAX, len(pred)) + 1)
                    cum = np.cumsum(pct[[pos[g] for g in pred[: kk[-1]]]]) / kk
                    topk_rows.append(pd.DataFrame({"universe": uname, "reference": rname, "TF": tf, "method": method,
                                                   "k": kk, "mean_response_percentile": cum}))
    return (pd.concat(run_rows, ignore_index=True), pd.concat(topk_rows, ignore_index=True),
            pd.concat(band_rows, ignore_index=True), pd.concat(hit_rows, ignore_index=True))


def method_style(method: str) -> dict:
    chrono = method.startswith("ChronoPack")
    return dict(color=b.METHOD_COLORS[method], lw=2.0 if chrono else 1.2, zorder=4 if chrono else 3, label=LEGEND.get(method, method))


def plot_topk_curves(topk: pd.DataFrame, universe: str) -> None:
    """Macro top-k response curves: mean over TFs; solid while every scored TF of the method has >= k predictions,
    dashed beyond (fewer TFs)."""
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), sharey=True)
    sub = topk[topk["universe"] == universe]
    for ax, ref in zip(axes, REFERENCES):
        s = sub[sub["reference"] == ref]
        for method in METHODS:
            d = s[s["method"] == method]
            if d.empty:
                continue
            n_tfs = d["TF"].nunique()
            g = d.groupby("k").agg(y=("mean_response_percentile", "mean"), n=("TF", "nunique")).reset_index()
            full = g[g["n"] == n_tfs]
            st = method_style(method)
            ax.plot(full["k"], full["y"], **st)
            tail = g[g["k"] >= full["k"].max()]
            if len(tail) > 1:
                ax.plot(tail["k"], tail["y"], ls=(0, (3, 2)), **{**st, "label": None, "lw": st["lw"] * 0.8})
        ax.axhline(0.5, color="#888888", lw=0.8, ls=":", zorder=1)
        ax.set_xscale("log")
        ax.set_xlim(1, CURVE_KMAX)
        ax.set_title(REFERENCES[ref], fontsize=10)
        ax.set_xlabel("top-k predicted targets per TF", fontsize=9.5)
        ax.tick_params(labelsize=8.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("mean response percentile of top-k\n(0.5 = random)", fontsize=9.5)
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.99, 0.5), frameon=False, fontsize=8.5)
    fig.suptitle(f"Perturbation response of predicted targets, averaged over TFs ({UNIVERSE_LABELS[universe]})", fontsize=10.5)
    fig.text(0.5, -0.02, "Response = -log10 P of the TF perturbation vs control (direction-agnostic). Dashed: k beyond "
             "some TFs' predicted-target count (mean over fewer TFs).", ha="center", fontsize=7.5, color="#555555")
    fig.tight_layout()
    stem = FIG_OUT / f"perturbseq_topk_response_curves_{universe}"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_gsea_curves(run: pd.DataFrame, band: pd.DataFrame, hits: pd.DataFrame, m: pd.DataFrame, counts: pd.DataFrame,
                     universe: str, ref: str) -> None:
    """Small multiples, one panel per TF: GSEA running enrichment of each method's top-K_MAIN set, over the 5-95%
    band (and median) of random sets of the same size; below each panel a hit rug (one row per method: where its
    top-K_MAIN targets fall in the response ranking), as in GSEA enrichment plots."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    s = run[(run["universe"] == universe) & (run["reference"] == ref)]
    bd = band[(band["universe"] == universe) & (band["reference"] == ref)]
    ht = hits[(hits["universe"] == universe) & (hits["reference"] == ref)]
    tfs = sorted(s["TF"].unique())
    ncol = 4
    nrow = int(np.ceil(len(tfs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(11.5, 3.3 * nrow), squeeze=False)
    c = counts[(counts["universe"] == universe) & (counts["reference"] == ref)].set_index("TF")
    zz = m[(m["universe"] == universe) & (m["reference"] == ref) & (m["K"] == K_MAIN)].set_index(["TF", "method"])
    for ax, tf in zip(axes.flat, tfs):
        d = s[s["TF"] == tf]
        bt = bd[bd["TF"] == tf]
        ax.fill_between(bt["rank_fraction"], bt["q05"], bt["q95"], color="#d9d9d9", lw=0, zorder=1)
        ax.plot(bt["rank_fraction"], bt["q50"], color="#999999", lw=0.9, ls=(0, (3, 2)), zorder=2)
        best = None
        for method in METHODS:
            dm = d[d["method"] == method]
            if dm.empty:
                continue
            ax.plot(dm["rank_fraction"], dm["running_ES"], **method_style(method))
            z = zz.loc[(tf, method), "ES_z"] if (tf, method) in zz.index else np.nan
            if np.isfinite(z) and (best is None or z > best[1]):
                best = (method, z)
        ax.axhline(0, color="#888888", lw=0.7)
        n_pos = int(c.loc[tf, "n_positives"]) if tf in c.index else 0
        ax.set_title(f"{tf}  ({n_pos} DEGs at BH FDR < 0.05)", fontsize=9)
        if best:
            ax.text(0.98, 0.95, f"best z: {LEGEND.get(best[0], best[0])} {best[1]:.1f}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=6.8, color="#444444")
        ax.tick_params(labelsize=7.5, labelbottom=False)
        ax.set_xlim(0, 1)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        rug = make_axes_locatable(ax).append_axes("bottom", size="42%", pad=0.03, sharex=ax)
        for i, method in enumerate(METHODS):
            hm = ht[(ht["TF"] == tf) & (ht["method"] == method)]
            if len(hm):
                rug.vlines(hm["rank_fraction"], i + 0.1, i + 0.9, color=b.METHOD_COLORS[method], lw=0.45)
        rug.set_ylim(len(METHODS), 0)
        rug.set_yticks([])
        rug.tick_params(labelsize=7.5)
        for sp in ("top", "right", "left"):
            rug.spines[sp].set_visible(False)
    for ax in axes.flat[len(tfs):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("running ES", fontsize=8.5)
    fig.supxlabel("genes ranked by perturbation response, strongest first (fraction of universe); rug: each method's "
                  f"top {K_MAIN} targets, one row per method in legend order", fontsize=9, y=0.02)
    handles = [plt.Line2D([], [], **{k: v for k, v in method_style(mm).items() if k != "zorder"}) for mm in METHODS]
    handles += [plt.Rectangle((0, 0), 1, 1, color="#d9d9d9", label=f"random {K_MAIN}-gene sets, 5-95%"),
                plt.Line2D([], [], color="#999999", lw=0.9, ls=(0, (3, 2)), label="random sets, median")]
    fig.legend(handles, [h.get_label() for h in handles], loc="lower center", ncol=6, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.06 if nrow > 2 else -0.1))
    fig.suptitle(f"GSEA running enrichment of each method's top {K_MAIN} targets: {REFERENCES[ref]} "
                 f"({UNIVERSE_LABELS[universe]})", fontsize=10.5)
    fig.tight_layout()
    stem = FIG_OUT / f"perturbseq_gsea_curves_{ref}_{universe}"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_curves(refs, edges, universes, m: pd.DataFrame, counts: pd.DataFrame) -> None:
    run, topk, band, hits = response_curves(refs, edges, universes)
    # the running-sum extreme must equal the ES in the metrics table
    es = run.loc[run.groupby(["universe", "reference", "TF", "method"])["running_ES"].apply(lambda v: v.abs().idxmax())]
    chk = es.merge(m[m["K"] == K_MAIN], on=["universe", "reference", "TF", "method"])
    assert len(chk) == len(m[m["K"] == K_MAIN]) and np.allclose(chk["running_ES"], chk["ES"]), "curve ES != table ES"
    run.to_csv(TAB_OUT / "perturbseq_gsea_running_curves.csv", index=False)
    topk.to_csv(TAB_OUT / "perturbseq_topk_response_curves.csv", index=False)
    band.to_csv(TAB_OUT / "perturbseq_gsea_random_set_band.csv", index=False)
    hits.to_csv(TAB_OUT / "perturbseq_gsea_hits.csv", index=False)
    for universe in universes:
        plot_topk_curves(topk, universe)
        for ref in REFERENCES:
            plot_gsea_curves(run, band, hits, m, counts, universe, ref)


ROC_GRID = np.linspace(0, 1, 501)
IN_CONTEXT_CHRONO = {"Replogle_K562": "ChronoPack_K562", "Yao_THP1": "ChronoPack_THP1"}
N_BOOT = 2000


def roc(scores: np.ndarray, labels: np.ndarray, n_pos: int, n_neg: int) -> tuple[np.ndarray, np.ndarray, float]:
    """ROC of a ranking over a universe of n_pos + n_neg items of which only `scores` / `labels` are ranked; the
    unranked items come last in random order (expected curve: a straight segment to (1, 1)), as in evaluate_one."""
    order = np.argsort(-scores, kind="mergesort")
    s, y = scores[order], labels[order].astype(bool)
    last = np.r_[np.flatnonzero(np.diff(s) != 0), len(s) - 1]  # tied scores are one step
    tp, fp = np.cumsum(y)[last], np.cumsum(~y)[last]
    tpr = np.r_[0.0, tp / n_pos, 1.0]
    fpr = np.r_[0.0, fp / n_neg, 1.0]
    return fpr, tpr, float(np.trapz(tpr, fpr))


def roc_curves(refs, edges_df, universes) -> tuple[pd.DataFrame, pd.DataFrame]:
    """BH FDR < 0.05 DEGs as labels (scMINER Fig. 4c / 4g style): per TF (>= MIN_POSITIVES positives) and pooled over all
    of a method's TFs in the reference (global network ROC; the pooled ranking uses the methods' edge weights)."""
    aucs, pts = [], []
    for uname, valid in universes.items():
        for rname, ref in refs.items():
            ref = ref[ref["target"].isin(valid)]
            pos_pairs = set(zip(ref.loc[ref["q_bh"] < FDR, "TF"], ref.loc[ref["q_bh"] < FDR, "target"]))
            measured = {tf: set(g["target"]) for tf, g in ref.groupby("TF")}
            e_all = edges_df[uname]
            for method in METHODS:
                e = e_all[(e_all["method"] == method) & e_all["TF"].isin(measured)]
                pooled_s, pooled_y, pooled_pos, pooled_n, tfs_used = [], [], 0, 0, []
                for tf, genes in measured.items():
                    et = e[(e["TF"] == tf) & e["target"].isin(genes)]
                    if et.empty:
                        continue
                    n_pos = sum((tf, g) in pos_pairs for g in genes)
                    y = np.array([(tf, g) in pos_pairs for g in et["target"]])
                    pooled_s.append(et["EdgeWeight"].to_numpy(float)); pooled_y.append(y)
                    pooled_pos += n_pos; pooled_n += len(genes); tfs_used.append(tf)
                    if n_pos >= MIN_POSITIVES:
                        # per TF: rank order only (ties in EdgeWeight kept as ties)
                        fpr, tpr, auc = roc(et["EdgeWeight"].to_numpy(float), y, n_pos, len(genes) - n_pos)
                        aucs.append(dict(universe=uname, reference=rname, scope="per_TF", TF=tf, method=method, AUROC=auc,
                                         n_positives=n_pos, n_universe=len(genes), n_ranked=len(et)))
                        pts.append(pd.DataFrame({"universe": uname, "reference": rname, "scope": "per_TF", "TF": tf,
                                                 "method": method, "FPR": ROC_GRID, "TPR": np.interp(ROC_GRID, fpr, tpr)}))
                if pooled_pos:
                    fpr, tpr, auc = roc(np.concatenate(pooled_s), np.concatenate(pooled_y), pooled_pos, pooled_n - pooled_pos)
                    aucs.append(dict(universe=uname, reference=rname, scope="global", TF=",".join(sorted(tfs_used)),
                                     method=method, AUROC=auc, n_positives=pooled_pos, n_universe=pooled_n,
                                     n_ranked=int(sum(len(v) for v in pooled_s))))
                    pts.append(pd.DataFrame({"universe": uname, "reference": rname, "scope": "global", "TF": "ALL",
                                             "method": method, "FPR": ROC_GRID, "TPR": np.interp(ROC_GRID, fpr, tpr)}))
    return pd.DataFrame(aucs), pd.concat(pts, ignore_index=True)


def check_roc(aucs: pd.DataFrame, refs, edges_df, universes) -> None:
    """Per-TF AUROC = Mann-Whitney U with the unranked genes tied at the bottom."""
    for _, row in aucs[aucs["scope"] == "per_TF"].iloc[::7].iterrows():
        r = refs[row["reference"]]
        r = r[(r["TF"] == row["TF"]) & r["target"].isin(universes[row["universe"]])]
        e = edges_df[row["universe"]]
        w = e[(e["method"] == row["method"]) & (e["TF"] == row["TF"])].set_index("target")["EdgeWeight"]
        sc = r["target"].map(w).fillna(-np.inf).to_numpy()
        sc[~np.isfinite(sc)] = np.nanmin(sc[np.isfinite(sc)]) - 1
        y = (r["q_bh"] < FDR).to_numpy()
        u = mannwhitneyu(sc[y], sc[~y], alternative="two-sided").statistic / (y.sum() * (~y).sum())
        assert abs(u - row["AUROC"]) < 1e-9, (row.to_dict(), u)


COVERAGE_NOTE = ("Full universe: DEGs outside the genes a method ranks fall in the tail after all of its ranked genes, so "
                 "AUROC < 0.5 here mostly reflects gene coverage (ChronoPack ranks ~9.5-10k genes, BEELINE ~2.3k); "
                 "see the shared-universe version.")


def plot_roc(aucs: pd.DataFrame, pts: pd.DataFrame, universe: str) -> None:
    pu, au = pts[pts["universe"] == universe], aucs[aucs["universe"] == universe]
    # global ROC, one panel per reference
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.3))
    for ax, ref in zip(axes, REFERENCES):
        a = au[(au["reference"] == ref) & (au["scope"] == "global")].set_index("method")
        for method in METHODS:
            d = pu[(pu["reference"] == ref) & (pu["scope"] == "global") & (pu["method"] == method)]
            if d.empty:
                continue
            st = method_style(method)
            st["label"] = f"{st['label']}  {a.loc[method, 'AUROC']:.2f} ({len(a.loc[method, 'TF'].split(','))} TFs)"
            ax.plot(d["FPR"], d["TPR"], **st)
        ax.plot([0, 1], [0, 1], color="#888888", lw=0.8, ls=":")
        n_pos = int(a["n_positives"].max()) if len(a) else 0
        ax.set_title(f"{REFERENCES[ref]}\n(up to {n_pos} TF-gene pairs with BH FDR < 0.05)", fontsize=9.5)
        ax.set_xlabel("false positive rate", fontsize=9)
        ax.legend(loc="lower right", fontsize=6.8, frameon=False, title="method  AUROC (TFs pooled)", title_fontsize=7)
        ax.set_aspect("equal")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("true positive rate", fontsize=9)
    fig.suptitle(f"Global ROC: all TF-gene pairs of each method's TFs, BH FDR < 0.05 DEGs as positives ({UNIVERSE_LABELS[universe]})", fontsize=10)
    fig.text(0.5, -0.01, "Pairs a method does not rank are placed last in random order (straight segment to (1, 1)). TF sets "
             "differ between methods (missing regulators)." + ("\n" + COVERAGE_NOTE if universe == "full" else ""),
             ha="center", va="top", fontsize=7.5, color="#555555")
    fig.tight_layout()
    stem = FIG_OUT / f"perturbseq_roc_global_{universe}"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # per TF, TFs with >= MIN_POSITIVES positives
    per = au[au["scope"] == "per_TF"]
    keys = per[["reference", "TF", "n_positives"]].drop_duplicates(["reference", "TF"]).sort_values(["reference", "n_positives"], ascending=[True, False])
    if keys.empty:
        return
    ncol = min(4, len(keys))
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.9 * ncol, 2.9 * nrow + 1.9), squeeze=False)
    for ax, (_, k) in zip(axes.flat, keys.iterrows()):
        a = per[(per["reference"] == k["reference"]) & (per["TF"] == k["TF"])].set_index("method")
        lines = []
        for method in METHODS:
            d = pu[(pu["reference"] == k["reference"]) & (pu["scope"] == "per_TF") & (pu["TF"] == k["TF"]) & (pu["method"] == method)]
            if d.empty:
                continue
            ax.plot(d["FPR"], d["TPR"], **method_style(method))
            lines.append((a.loc[method, "AUROC"], method))
        ax.plot([0, 1], [0, 1], color="#888888", lw=0.8, ls=":")
        for j, (v, method) in enumerate(sorted(lines, reverse=True)):
            ax.text(0.98, 0.03 + 0.065 * j, f"{SHORT.get(method, method)} {v:.2f}", transform=ax.transAxes, ha="right",
                    va="bottom", fontsize=5.8, color=b.METHOD_COLORS[method])
        ctx = "K562" if k["reference"] == "Replogle_K562" else "THP1"
        ax.set_title(f"{k['TF']} | {ctx} perturb ({int(k['n_positives'])} DEGs)", fontsize=8.5)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    for ax in axes.flat[len(keys):]:
        ax.axis("off")
    handles = [plt.Line2D([], [], **{k: v for k, v in method_style(mm).items() if k != "zorder"}) for mm in METHODS]
    finish_grid(fig, handles, nrow, "false positive rate", "true positive rate",
                f"Per-TF ROC, BH FDR < 0.05 DEGs as positives ({UNIVERSE_LABELS[universe]}; TFs with >= {MIN_POSITIVES} DEGs; "
                "AUROC per panel, Chrono = ChronoPack, K562 / THP1 by colour)",
                "Genes a method does not rank are placed last in random order (straight segment to (1, 1))."
                + ("\n" + COVERAGE_NOTE if universe == "full" else ""))
    stem = FIG_OUT / f"perturbseq_roc_per_tf_{universe}"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def finish_grid(fig, handles, nrow: int, xlabel: str, ylabel: str, title: str, note: str) -> None:
    """Small-multiple layout with fixed margins (in inches) for the title, shared axis labels, legend and footnote."""
    w, h = fig.get_size_inches()
    n_lines = note.count("\n") + 1
    top, bottom = 0.55, 1.05 + 0.14 * n_lines  # footnote lines + legend (2 rows) + x label
    fig.tight_layout(rect=(0.35 / w, bottom / h, 1, 1 - top / h))
    fig.suptitle(title, fontsize=9.5, y=1 - 0.12 / h, va="top")
    fig.supxlabel(xlabel, fontsize=9, y=(bottom - 0.42) / h, va="bottom")
    fig.supylabel(ylabel, fontsize=9, x=0.05 / w, ha="left")
    fig.legend(handles, [hd.get_label() for hd in handles], loc="upper center", ncol=5, frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, (bottom - 0.47) / h))
    fig.text(0.5, 0.03 / h, note, ha="center", va="bottom", fontsize=7.2, color="#555555")


SWEEP_N = [10, 25, 50, 100, 200, 500]
MIN_SOFT = 2.0  # soft-label ROC drawn when the universe holds >= this much responder mass (sum of 1 - lfdr)


def roc_soft(scores: np.ndarray, soft: np.ndarray, tot_pos: float, tot_neg: float) -> tuple[np.ndarray, np.ndarray, float]:
    """ROC with soft labels: each gene counts as `soft` of a positive and 1 - `soft` of a negative; TPR / FPR are the
    cumulative positive / negative mass down the ranking (tied scores one step), unranked genes last in random order.
    With 0/1 labels this is the ordinary ROC (`roc`)."""
    order = np.argsort(-scores, kind="mergesort")
    s, w = scores[order], soft[order]
    last = np.r_[np.flatnonzero(np.diff(s) != 0), len(s) - 1]
    tpr = np.r_[0.0, np.cumsum(w)[last] / tot_pos, 1.0]
    fpr = np.r_[0.0, np.cumsum(1 - w)[last] / tot_neg, 1.0]
    return fpr, tpr, float(np.trapz(tpr, fpr))


def weighted_mw(scores: np.ndarray, soft: np.ndarray) -> float:
    """Independent check of roc_soft: P(responder mass scores above non-responder mass), ties (incl. a gene with
    itself) count 1/2; scores over the whole universe (unranked genes share the lowest score)."""
    df = pd.DataFrame({"s": scores, "p": soft, "n": 1 - soft}).groupby("s")[["p", "n"]].sum().sort_index()
    below = np.r_[0.0, np.cumsum(df["n"].to_numpy())[:-1]]
    return float((df["p"].to_numpy() * (below + 0.5 * df["n"].to_numpy())).sum() / (df["p"].sum() * df["n"].sum()))


def soft_and_sweep(refs, edges_df, universes) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Per (universe, reference, TF, method): soft-label AUROC and curve (1 - lfdr as labels), and the threshold sweep
    (binary AUROC with the top N responders by P as positives)."""
    soft_rows, soft_pts, sweep_rows = [], [], []
    for uname, valid in universes.items():
        for rname, ref in refs.items():
            e_all = edges_df[uname]
            for tf, r in ref.groupby("TF", sort=False):
                r = r[r["target"].isin(valid)]
                r = r.sort_values(["score", "target"], ascending=[False, True], kind="mergesort").reset_index(drop=True)
                genes = r["target"].to_numpy()
                soft = r["soft"].to_numpy()
                tot_pos, tot_neg = soft.sum(), (1 - soft).sum()
                for method in METHODS:
                    w = e_all[(e_all["method"] == method) & (e_all["TF"] == tf)].set_index("target")["EdgeWeight"]
                    w = w[w.index.isin(set(genes))]
                    if w.empty:
                        continue
                    ranked = np.isin(genes, w.index)
                    sc = r.loc[ranked, "target"].map(w).to_numpy(float)
                    if tot_pos >= MIN_SOFT:
                        fpr, tpr, auc = roc_soft(sc, soft[ranked], tot_pos, tot_neg)
                        soft_rows.append(dict(universe=uname, reference=rname, TF=tf, method=method, soft_AUROC=auc,
                                              responder_mass=tot_pos, pi0=float(r["pi0"].iloc[0]), n_universe=len(genes),
                                              n_ranked=int(ranked.sum())))
                        soft_pts.append(pd.DataFrame({"universe": uname, "reference": rname, "TF": tf, "method": method,
                                                      "FPR": ROC_GRID, "TPR": np.interp(ROC_GRID, fpr, tpr)}))
                    for n_top in SWEEP_N:
                        if n_top >= len(genes):
                            continue
                        y = np.zeros(len(genes), dtype=bool)
                        y[:n_top] = True  # r is ordered by response
                        _, _, auc = roc(sc, y[ranked], n_top, len(genes) - n_top)
                        sweep_rows.append(dict(universe=uname, reference=rname, TF=tf, method=method, N=n_top, AUROC=auc,
                                               n_bh05=int((r["q_bh"] < FDR).sum())))
    return pd.DataFrame(soft_rows), pd.concat(soft_pts, ignore_index=True), pd.DataFrame(sweep_rows)


def check_soft(sa: pd.DataFrame, refs, edges_df, universes) -> None:
    for _, row in sa.iloc[::9].iterrows():
        r = refs[row["reference"]]
        r = r[(r["TF"] == row["TF"]) & r["target"].isin(universes[row["universe"]])]
        e = edges_df[row["universe"]]
        w = e[(e["method"] == row["method"]) & (e["TF"] == row["TF"])].set_index("target")["EdgeWeight"]
        sc = r["target"].map(w).to_numpy(float)
        sc[~np.isfinite(sc)] = np.nanmin(sc) - 1
        u = weighted_mw(sc, r["soft"].to_numpy())
        assert abs(u - row["soft_AUROC"]) < 1e-9, (row.to_dict(), u)


def plot_soft_roc(sa: pd.DataFrame, sp: pd.DataFrame, counts: pd.DataFrame, universe: str) -> None:
    sa, sp = sa[sa["universe"] == universe], sp[sp["universe"] == universe]
    keys = (sa[["reference", "TF", "responder_mass", "pi0"]].drop_duplicates(["reference", "TF"])
            .sort_values(["reference", "responder_mass"], ascending=[True, False]))
    if keys.empty:
        return
    c = counts[counts["universe"] == universe].set_index(["reference", "TF"])
    ncol = min(4, len(keys))
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.9 * ncol, 2.95 * nrow + 2.2), squeeze=False)
    for ax, (_, k) in zip(axes.flat, keys.iterrows()):
        a = sa[(sa["reference"] == k["reference"]) & (sa["TF"] == k["TF"])].set_index("method")
        lines = []
        for method in METHODS:
            d = sp[(sp["reference"] == k["reference"]) & (sp["TF"] == k["TF"]) & (sp["method"] == method)]
            if d.empty:
                continue
            ax.plot(d["FPR"], d["TPR"], **method_style(method))
            lines.append((a.loc[method, "soft_AUROC"], method))
        ax.plot([0, 1], [0, 1], color="#888888", lw=0.8, ls=":")
        for j, (v, method) in enumerate(sorted(lines, reverse=True)):
            ax.text(0.98, 0.03 + 0.065 * j, f"{SHORT.get(method, method)} {v:.2f}", transform=ax.transAxes, ha="right",
                    va="bottom", fontsize=5.8, color=b.METHOD_COLORS[method])
        ctx = "K562" if k["reference"] == "Replogle_K562" else "THP1"
        n_bh = int(c.loc[(k["reference"], k["TF"]), "n_positives"]) if (k["reference"], k["TF"]) in c.index else 0
        ax.set_title(f"{k['TF']} | {ctx} perturb\nresponder mass {k['responder_mass']:.0f} (BH 0.05: {n_bh})", fontsize=8)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=7)
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
    for ax in axes.flat[len(keys):]:
        ax.axis("off")
    handles = [plt.Line2D([], [], **{k: v for k, v in method_style(mm).items() if k != "zorder"}) for mm in METHODS]
    xl, yl = "false positive rate (cumulative non-responder mass, lfdr)", "true positive rate (cumulative responder mass, 1 - lfdr)"
    omitted = sorted(set(map(tuple, counts[(counts["universe"] == universe) & (counts["soft_positives"] < MIN_SOFT)][["reference", "TF"]].to_numpy())))
    om = ", ".join(f"{t} ({'K562' if r_ == 'Replogle_K562' else 'THP1'})" for r_, t in omitted)
    note = ("Local FDR from a beta-uniform mixture fit to each TF's P values (Pounds & Morris 2003). Responder mass = "
             "expected number of responding genes in the universe. Genes a method does not rank are last in random order."
             + (f"\nOmitted (responder mass < {MIN_SOFT:g}, no detectable response): {om}." if om else "")
             + ("\n" + COVERAGE_NOTE if universe == "full" else ""))
    finish_grid(fig, handles, nrow, xl, yl, "Per-TF soft-label ROC, threshold-free: each gene weighted by its probability of "
                f"responding, 1 - local FDR ({UNIVERSE_LABELS[universe]})", note)
    stem = FIG_OUT / f"perturbseq_roc_soft_per_tf_{universe}"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_sweep(sw: pd.DataFrame, universe: str) -> None:
    sw = sw[sw["universe"] == universe]
    keys = sw[["reference", "TF"]].drop_duplicates().sort_values(["reference", "TF"])
    ncol = 6
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.25 * ncol, 2.0 * nrow + 0.7), sharex=True, sharey=True, squeeze=False)
    for ax, (_, k) in zip(axes.flat, keys.iterrows()):
        d = sw[(sw["reference"] == k["reference"]) & (sw["TF"] == k["TF"])]
        for method in METHODS:
            dm = d[d["method"] == method]
            if len(dm):
                ax.plot(dm["N"], dm["AUROC"], marker="o", ms=2.2, **method_style(method))
        ax.axhline(0.5, color="#888888", lw=0.7, ls=":")
        n_bh = int(d["n_bh05"].iloc[0])
        if n_bh >= SWEEP_N[0]:
            ax.axvline(n_bh, color="#bbbbbb", lw=0.8, ls="--")
        ax.set_xscale("log")
        ctx = "K562" if k["reference"] == "Replogle_K562" else "THP1"
        ax.set_title(f"{k['TF']} | {ctx}", fontsize=8)
        ax.tick_params(labelsize=6.5)
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
    for ax in axes.flat[len(keys):]:
        ax.axis("off")
    handles = [plt.Line2D([], [], **{k: v for k, v in method_style(mm).items() if k != "zorder"}) for mm in METHODS]
    handles.append(plt.Line2D([], [], color="#bbbbbb", lw=0.8, ls="--", label="BH FDR < 0.05 count (all measured genes)"))
    fig.legend(handles, [h.get_label() for h in handles], loc="lower center", ncol=5, frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, -0.06))
    fig.supxlabel("N most responsive genes (by P) used as positives", fontsize=9, y=0.02)
    fig.supylabel("AUROC (0.5 = random)", fontsize=9)
    fig.suptitle(f"Threshold sweep: per-TF AUROC with the top N responders as positives ({UNIVERSE_LABELS[universe]})", fontsize=10)
    fig.tight_layout()
    stem = FIG_OUT / f"perturbseq_roc_threshold_sweep_{universe}"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def summary_stats(m: pd.DataFrame) -> pd.DataFrame:
    """scMINER Fig. 4d style: per method, median over TFs with a TF-bootstrap 95% CI, and a paired Wilcoxon test
    (two-sided, over the TFs both have) against the in-context ChronoPack network."""
    from scipy.stats import wilcoxon

    km = m[m["K"] == K_MAIN].assign(neglog10_gsea_p=lambda d: -np.log10(d["perm_p"]))
    rows = []
    rng = np.random.default_rng(SEED)
    for (uname, rname), g in km.groupby(["universe", "reference"]):
        ref_m = IN_CONTEXT_CHRONO[rname]
        for metric in ("response_AUROC", "neglog10_gsea_p"):
            base = g[g["method"] == ref_m].set_index("TF")[metric]
            for method in METHODS:
                v = g[g["method"] == method].set_index("TF")[metric]
                if v.empty:
                    continue
                boots = np.median(rng.choice(v.to_numpy(), size=(N_BOOT, len(v)), replace=True), axis=1)
                lo, hi = np.percentile(boots, [2.5, 97.5])
                shared = v.index.intersection(base.index)
                p = np.nan
                if method != ref_m and len(shared) >= 5 and (v[shared] != base[shared]).any():
                    p = float(wilcoxon(v[shared], base[shared], alternative="two-sided").pvalue)
                rows.append(dict(universe=uname, reference=rname, metric=metric, method=method, n_TFs=len(v),
                                 median=float(v.median()), ci_low=lo, ci_high=hi, compared_to=ref_m,
                                 n_paired=len(shared), wilcoxon_p=p))
    return pd.DataFrame(rows)


def plot_summary(km: pd.DataFrame, st: pd.DataFrame, universe: str) -> None:
    km = km[(km["universe"] == universe) & (km["K"] == K_MAIN)].assign(neglog10_gsea_p=lambda d: -np.log10(d["perm_p"]))
    metrics = {"response_AUROC": ("response AUROC, top-100 targets (0.5 = random)", 0.5),
               "neglog10_gsea_p": ("-log10 GSEA permutation P (cap 3.0)", -np.log10(0.05))}
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2))
    rng = np.random.default_rng(1)
    for i, ref in enumerate(REFERENCES):
        for j, (metric, (ylabel, line)) in enumerate(metrics.items()):
            ax = axes[i, j]
            s = st[(st["universe"] == universe) & (st["reference"] == ref) & (st["metric"] == metric)].set_index("method")
            ms = [mm for mm in METHODS if mm in s.index]
            x = np.arange(len(ms))
            ax.bar(x, s.loc[ms, "median"], color=[b.METHOD_COLORS[mm] for mm in ms], width=0.65, alpha=0.85, zorder=2)
            ax.errorbar(x, s.loc[ms, "median"], yerr=[s.loc[ms, "median"] - s.loc[ms, "ci_low"], s.loc[ms, "ci_high"] - s.loc[ms, "median"]],
                        fmt="none", ecolor="#333333", elinewidth=0.9, capsize=2.5, zorder=3)
            for xi, mm in zip(x, ms):
                v = km[(km["reference"] == ref) & (km["method"] == mm)][metric].to_numpy()
                ax.scatter(xi + rng.uniform(-0.18, 0.18, len(v)), v, s=7, color="#222222", alpha=0.55, lw=0, zorder=4)
                p = s.loc[mm, "wilcoxon_p"]
                if np.isfinite(p):
                    ax.text(xi, 1.0, f"P={p:.2g}", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.3,
                            color="#444444")
            ax.axhline(line, color="#888888", lw=0.8, ls=":", zorder=1)
            if metric == "response_AUROC":
                lo = min(0.3, float(km[(km["reference"] == ref)][metric].min()) - 0.02)
                ax.set_ylim(lo, max(0.82, float(km[(km["reference"] == ref)][metric].max()) + 0.03))
            else:
                ax.set_ylim(0, 3.2)
            ax.set_xticks(x)
            ax.set_xticklabels([LEGEND.get(mm, mm).replace("ChronoPack | ", "Chrono\n") for mm in ms], rotation=45, ha="right", fontsize=7.5)
            ax.tick_params(axis="y", labelsize=7.5)
            ax.set_ylabel(ylabel, fontsize=8)
            ax.set_title(f"{REFERENCES[ref]}", fontsize=9, pad=14)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
    fig.suptitle(f"Per-TF perturbation benchmark summary ({UNIVERSE_LABELS[universe]})", fontsize=10.5)
    fig.text(0.5, -0.01, "Bars: median over TFs; error bars: 95% TF-bootstrap CI; dots: TFs. P: paired two-sided Wilcoxon vs "
             "the in-context ChronoPack (K562 for Replogle, THP1 for Yao) over shared TFs (>= 5). Dotted: random (AUROC) / P = 0.05.",
             ha="center", fontsize=7.3, color="#555555")
    fig.tight_layout()
    stem = FIG_OUT / f"perturbseq_summary_{universe}"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_benchmark_panels(refs, edges_df, universes, m: pd.DataFrame, counts: pd.DataFrame) -> None:
    aucs, pts = roc_curves(refs, edges_df, universes)
    check_roc(aucs, refs, edges_df, universes)
    sa, sp, sw = soft_and_sweep(refs, edges_df, universes)
    check_soft(sa, refs, edges_df, universes)
    sa.to_csv(TAB_OUT / "perturbseq_roc_soft_auc.csv", index=False)
    sp.to_csv(TAB_OUT / "perturbseq_roc_soft_curves.csv", index=False)
    sw.to_csv(TAB_OUT / "perturbseq_roc_threshold_sweep.csv", index=False)
    st = summary_stats(m)
    aucs.to_csv(TAB_OUT / "perturbseq_roc_auc.csv", index=False)
    pts.to_csv(TAB_OUT / "perturbseq_roc_curves.csv", index=False)
    st.to_csv(TAB_OUT / "perturbseq_summary_stats.csv", index=False)
    for universe in universes:
        plot_roc(aucs, pts, universe)
        plot_soft_roc(sa, sp, counts, universe)
        plot_sweep(sw, universe)
        plot_summary(m, st, universe)


def checks(m: pd.DataFrame, refs, edges, universes) -> None:
    # Mann-Whitney AUROC against scipy for every 25th row
    for _, row in m.iloc[::25].iterrows():
        r = refs[row["reference"]]
        r = r[(r["TF"] == row["TF"]) & r["target"].isin(universes[row["universe"]])]
        pred = [g for g in edges[row["universe"]][(row["method"], row["TF"])] if g in set(r["target"])][: row["K"]]
        a = r[r["target"].isin(pred)]["score"]
        o = r[~r["target"].isin(pred)]["score"]
        u = mannwhitneyu(a, o, alternative="two-sided").statistic / (len(a) * len(o))
        assert abs(u - row["response_AUROC"]) < 1e-9, (row.to_dict(), u)
    assert not set(m["method"]) & set(fig5_paths.PERTURB_METHODS)
    assert (m[m["universe"] == "shared"]["n_universe"] <= len(universes["shared"])).all()


def main() -> None:
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    TAB_OUT.mkdir(parents=True, exist_ok=True)
    refs = load_references()
    universes = {"shared": set(SHARED_GENES.read_text().split()), "full": b.load_valid_targets()}
    edges_df = {}
    for u, p in EDGES.items():
        e = pd.read_csv(p, usecols=["method", "TF", "target", "EdgeWeight"])
        edges_df[u] = e[e["method"].isin(METHODS)]
    edges = {u: load_edges(u) for u in EDGES}

    if "--curves-only" in sys.argv:  # redraw curves, ROC and summary panels from the existing metric tables
        m = pd.read_csv(TAB_OUT / "perturbseq_rank_metrics.csv")
        make_curves(refs, edges, universes, m, pd.read_csv(TAB_OUT / "perturbseq_reference_positives.csv"))
        make_benchmark_panels(refs, edges_df, universes, m, pd.read_csv(TAB_OUT / "perturbseq_reference_positives.csv"))
        print("PERTURBSEQ_CURVES_OK")
        return
    m = rank_metrics(refs, edges, universes)
    checks(m, refs, edges, universes)
    f, counts = fdr_metrics(refs, edges_df, universes)
    m.to_csv(TAB_OUT / "perturbseq_rank_metrics.csv", index=False)
    f.to_csv(TAB_OUT / "perturbseq_fdr05_metrics.csv", index=False)
    counts.to_csv(TAB_OUT / "perturbseq_reference_positives.csv", index=False)

    for universe in universes:
        for metric in ("response_AUROC", "ES_z"):
            plot_rank(m, universe, metric)
        plot_fdr(f, counts, universe)
    make_curves(refs, edges, universes, m, counts)
    make_benchmark_panels(refs, edges_df, universes, m, counts)

    km = m[m["K"] == K_MAIN]
    med = (km.groupby(["universe", "reference", "method"])
           .agg(n_TFs=("TF", "size"), median_response_AUROC=("response_AUROC", "median"), median_ES_z=("ES_z", "median"), median_NES=("NES", "median"),
                n_sig_q05=("perm_q", lambda q: int((q < 0.05).sum())))
           .reset_index())
    med.to_csv(TAB_OUT / "perturbseq_rank_medians_K100.csv", index=False)
    summary = {"references": {k: sorted(v["TF"].unique()) for k, v in refs.items()},
               "methods": METHODS, "K": K_MAIN, "n_perm": N_PERM, "universe_sizes": {k: len(v) for k, v in universes.items()},
               "rank_rows": int(len(m)), "fdr_rows": int(len(f))}
    (TAB_OUT / "perturbseq_reference_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    pd.set_option("display.width", 200)
    print(counts.pivot_table(index=["reference", "TF"], columns="universe", values="n_positives").to_string())
    print(med.to_string(index=False))
    print(f"Wrote {FIG_OUT} and {TAB_OUT}")
    print("PERTURBSEQ_REFERENCE_V2_OK")


if __name__ == "__main__":
    main()
