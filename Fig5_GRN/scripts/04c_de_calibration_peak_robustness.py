#!/usr/bin/env python
"""DE score calibration and peak-time robustness (side analysis; Fig5 outputs are not changed).

Part A -- calibration. Compares the distributions of ChronoPack's target scores (tradeSeq 1 - lfdr_bum)
and the DE-derived scores (04a min-p DE 1 - lfdr_bum, used by the 04 global ablation; peak DE), and re-runs
the TF-balanced global evaluation of 04 with variants that separate gene ORDER from score CALIBRATION:
  TF x tradeSeq 1-lfdr                      (full ChronoPack)
  TF x DE(min-p) 1-lfdr                     (the 04 "TF prior x DE" variant; reproduced)
  TF x DE(peak, fitted) 1-lfdr              (peak DE, bounded via ChronoPack's BUM calibration)
  TF x DE(min-p) order, tradeSeq calibration (quantile-matched: DE ranks get tradeSeq's sorted scores)
  TF x DE(peak) order, tradeSeq calibration
  TF x tradeSeq order, DE(min-p) calibration (tradeSeq ranks get DE's sorted scores)
Quantile matching keeps each gene ordering and replaces only the score values, so a difference between a
variant and its quantile-matched twin is calibration, and between twins with the same calibration is ordering.

Part B -- peak-time robustness. Peak time point per gene under three schemes:
  pcdf      pc_df.new.csv t_peak mapped from the aligned clustering axis (K562a: minutes, capped at 97 by the
            shared window; THP1: t_peak + 10/0.65), observed peak for genes outside pc_df (= 04b)
  fitted    peak of ChronoPack's fitted curve (screen fit sce_screen.1trajectory.nknots4.rds, all tested genes),
            native time axis from 0 (K562 0-110, THP1 0-165 min), |log2((mu+0.01)/(mu0+0.01))| maximal, no cap
  observed  time point with the largest |log2FC| of mean normalized expression vs baseline, all genes
Each is mapped to the nearest sampled time point; DE = Wilcoxon of that time point vs baseline (as in 04b).
Reports the time-axis check (pc_df conversion vs fitted peaks), the K562 97-min cap, genes absent from
pc_df, and per-TF EPR/AUPRCr (Fig5 evaluation) for each scheme.

Writes output/de_calibration_peak/{tables,figures}/.
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
import scanpy as sc

import chronopack as cp
import benchmark_grn_systems as b
import compare_regulator_prioritization as compare
import generate_system_broad_filtered_publication_figures_v2 as pub_v2


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name.replace("-", "_"), SCRIPTS_DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


g04 = _load("04_tf_balanced_global_recall")  # TF-balanced recall engine and settings
g04b = _load("04b_per_tf_de_peak")  # fit inputs, pc_df peak conversion

OUT = fig5_paths.OUTPUT_DIR / "de_calibration_peak"
TABLES, FIGS = OUT / "tables", OUT / "figures"
PROTEIN_CODING = fig5_paths.CHRONOPACK_INPUT_DIR / "protein_coding_genes.csv"
SCREEN_SCE = "sce_screen.1trajectory.nknots4.rds"
FIT_PSEUDOCOUNT = 0.01  # as cp.tl.signed_max_log2fc_from_start
SCHEMES = ["pcdf", "fitted", "observed"]
GLOBAL_VARIANTS = {
    "tradeseq": "TF x tradeSeq 1-lfdr (full)",
    "de_minp": "TF x DE(min-p) 1-lfdr (04 ablation)",
    "de_peak": "TF x DE(peak, fitted) 1-lfdr",
    "de_minp_qm": "DE(min-p) order, tradeSeq calibration",
    "de_peak_qm": "DE(peak) order, tradeSeq calibration",
    "tradeseq_qm_de": "tradeSeq order, DE(min-p) calibration",
}


def quantile_match(order_scores: pd.DataFrame, values: np.ndarray) -> pd.DataFrame:
    """Give the genes ranked by (score, tie) the sorted ``values`` (descending), keeping their tie statistic."""
    o = order_scores.sort_values(["score", "waldStat_1"], ascending=False, kind="mergesort").reset_index(drop=True)
    v = np.sort(np.asarray(values, float))[::-1][: len(o)]
    return o.assign(score=v)


def fitted_peaks(context: str, genes: list[str]) -> pd.DataFrame:
    _cell_type, minutes_per_tag, _conv, t_max = g04b.CONTEXT[context]
    fit_dir = fig5_paths.CHRONOPACK_INPUT_DIR / context / "tradeseq_fit"
    curves = cp.tl.predict_gene_dynamics(fit_dir, genes, n_points=int(t_max) + 1, time_max=t_max, sce_file=SCREEN_SCE)
    mu = np.expm1(curves.pivot(index="gene", columns="t", values="mu_log1p")).reindex(genes)
    lfc = np.log2((mu + FIT_PSEUDOCOUNT).div(mu.iloc[:, 0] + FIT_PSEUDOCOUNT, axis=0))
    later = lfc.loc[:, lfc.columns > 0]
    peak_min = later.abs().idxmax(axis=1).astype(float)
    peak_lfc = pd.Series([later.loc[g, m] for g, m in zip(later.index, peak_min)], index=later.index)
    tag = (np.round(peak_min / minutes_per_tag).astype(int) + 1).clip(lower=2, upper=int(round(t_max / minutes_per_tag)) + 1)
    return pd.DataFrame({"fitted_peak_min": peak_min, "fitted_peak_log2fc": peak_lfc, "fitted_peak_tag": tag})


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(g04.SEED)
    valid = b.load_valid_targets()
    gts = pub_v2.load_reference_ground_truths(valid)

    # -------------------------------------------------------------- peaks, DE and scores per context
    target_scores, diag_rows, dist_rows = {}, [], []
    peaks_all = {}
    for context in g04b.CONTEXT:
        a = g04b.load_fit(context)
        peaks = g04b.peak_tags(context, a)  # pc_df / observed scheme of 04b
        asso = pd.read_csv(fig5_paths.CHRONOPACK_INPUT_DIR / context / "asso_df.1trajectory.nknots4.csv", index_col=0)
        genes = [g for g in peaks.index if g in asso.index]
        peaks = peaks.loc[genes].join(fitted_peaks(context, genes))
        tags = list(a.obs["time"].cat.categories)
        sc.tl.rank_genes_groups(a, groupby="time", groups=tags[1:], reference=tags[0], method="wilcoxon", n_genes=a.n_vars)
        res = a.uns["rank_genes_groups"]
        p_tab = pd.DataFrame({int(t): pd.Series(res["pvals"][t], index=res["names"][t]) for t in tags[1:]}).reindex(genes)
        z_tab = pd.DataFrame({int(t): pd.Series(res["scores"][t], index=res["names"][t]) for t in tags[1:]}).reindex(genes).abs()
        scheme_tag = {"pcdf": peaks["peak_tag"], "fitted": peaks["fitted_peak_tag"], "observed": peaks["observed_peak_tag"]}
        for s, tg in scheme_tag.items():
            peaks[f"p_{s}"] = [p_tab.loc[g, t] for g, t in zip(genes, tg)]
            peaks[f"z_{s}"] = [z_tab.loc[g, t] for g, t in zip(genes, tg)]
        peaks_all[context] = peaks
        peaks.to_csv(TABLES / f"peaks_{context}.csv")

        # Time-axis and cap diagnostics.
        pcg = peaks["peak_source"] == "pc_df"
        cell_type, mpt, conv, t_max = g04b.CONTEXT[context]
        pc = pd.read_csv(g04b.PC_DF)
        pc = pc[pc["cell_type"] == cell_type].set_index("gene").reindex(peaks.index[pcg])
        pc_min = pc["t_peak"].map(conv)
        fitted_min = peaks.loc[pcg, "fitted_peak_min"]
        capped = pc["t_peak"] >= 97 - 1e-9 if context == "K562" else pd.Series(False, index=pc.index)
        unc = ~capped
        diag_rows.append({
            "context": context, "n_tested_genes": len(peaks), "n_in_pc_df": int(pcg.sum()), "n_absent_from_pc_df": int((~pcg).sum()),
            "pc_df_t_peak_min_context_min": float(pc_min.min()), "pc_df_t_peak_max_context_min": float(pc_min.max()),
            "n_capped_at_97": int(capped.sum()),
            "capped_fitted_peak_median_min": float(fitted_min[capped].median()) if capped.any() else np.nan,
            "pcdf_vs_fitted_minutes_median_diff_uncapped": float((pc_min[unc] - fitted_min[unc]).median()),
            "pcdf_vs_fitted_minutes_spearman_uncapped": float(pd.Series(pc_min[unc]).corr(fitted_min[unc], method="spearman")),
            "pcdf_vs_fitted_same_tag_uncapped": float((peaks.loc[pcg, "pc_df_peak_tag"][unc] == peaks.loc[pcg, "fitted_peak_tag"][unc]).mean()),
            "pcdf_vs_fitted_within1_tag_uncapped": float(((peaks.loc[pcg, "pc_df_peak_tag"][unc] - peaks.loc[pcg, "fitted_peak_tag"][unc]).abs() <= 1).mean()),
            "fitted_vs_observed_same_tag_all": float((peaks["fitted_peak_tag"] == peaks["observed_peak_tag"]).mean()),
            "fitted_vs_observed_within1_all": float(((peaks["fitted_peak_tag"] - peaks["observed_peak_tag"]).abs() <= 1).mean()),
            "sign_pcdf_vs_fitted": float((np.sign(pc["peak_log2fc"]) == np.sign(peaks.loc[pcg, "fitted_peak_log2fc"])).mean()),
            "fitted_peak_before_17min_in_pc_df_genes": int((fitted_min < 10 / 0.65 + 1).sum()) if context == "THP1" else np.nan,
        })

        # Scores: tradeSeq, DE(min-p) 1-lfdr (04a), DE(peak, each scheme) 1-lfdr and -log10 p.
        ts = pd.read_csv(fig5_paths.CHRONOPACK_OUTPUT_DIR / context / "target_scores.scores.csv")
        dm = pd.read_csv(g04.DE_DIR / context / "target_scores.scores.csv")
        sc_tab = {"tradeseq": ts.assign(score=1 - ts.lfdr_bum)[["gene", "score", "waldStat_1"]],
                  "de_minp": dm.assign(score=1 - dm.lfdr_bum)[["gene", "score", "waldStat_1"]]}
        for s in SCHEMES:
            assoc = TABLES / f"peak_de_association_{context}_{s}.csv"
            pd.DataFrame({"pvalue_1": peaks[f"p_{s}"], "waldStat_1": peaks[f"z_{s}"]}, index=peaks.index.rename("gene")).to_csv(assoc)
            lf = cp.tl.grn.score_targets(assoc, valid_targets=PROTEIN_CODING).scores
            sc_tab[f"de_peak_{s}"] = lf.assign(score=1 - lf.lfdr_bum)[["gene", "score", "waldStat_1"]]
            sc_tab[f"de_peak_{s}_logp"] = pd.DataFrame({"gene": peaks.index, "score": -np.log10(peaks[f"p_{s}"].clip(lower=1e-300)).to_numpy(),
                                                        "waldStat_1": peaks[f"z_{s}"].to_numpy()})
        for k in sc_tab:
            t = sc_tab[k]
            sc_tab[k] = t[t["gene"].astype(str).isin(valid) & np.isfinite(t["score"]) & (t["score"] > 0)].drop_duplicates("gene").reset_index(drop=True)
        common = set.intersection(*[set(t["gene"]) for t in sc_tab.values()])
        for k in sc_tab:
            sc_tab[k] = sc_tab[k][sc_tab[k]["gene"].isin(common)].reset_index(drop=True)
        sc_tab["de_minp_qm"] = quantile_match(sc_tab["de_minp"], sc_tab["tradeseq"]["score"])
        sc_tab["de_peak_qm"] = quantile_match(sc_tab["de_peak_fitted"], sc_tab["tradeseq"]["score"])
        sc_tab["tradeseq_qm_de"] = quantile_match(sc_tab["tradeseq"], sc_tab["de_minp"]["score"])
        sc_tab["de_peak"] = sc_tab["de_peak_fitted"]
        target_scores[context] = sc_tab
        for k in ["tradeseq", "de_minp", "de_peak_pcdf", "de_peak_fitted", "de_peak_observed"]:
            s = np.sort(sc_tab[k]["score"].to_numpy())[::-1]
            dist_rows.append({"context": context, "score": k, "n": len(s), "max": s[0], "median": float(np.median(s)),
                              **{f">{c}": int((s > c).sum()) for c in (0.5, 0.9, 0.99)},
                              **{f"at_rank_{r}": float(s[r - 1]) for r in (10, 100, 300, 1000, 3000)},
                              "frac_mass_top100": float(s[:100].sum() / s.sum())})
        print(f"{context}: {len(common)} genes scored by every variant", flush=True)

    diag = pd.DataFrame(diag_rows)
    dist = pd.DataFrame(dist_rows)
    diag.to_csv(TABLES / "peak_time_diagnostics.csv", index=False)
    dist.to_csv(TABLES / "score_distributions.csv", index=False)
    pd.set_option("display.width", 250)
    print(diag.T.to_string())
    print(dist.round(4).to_string(index=False), flush=True)

    # -------------------------------------------------------------- per-TF metrics per peak scheme
    pertf = []
    for context, sc_tab in target_scores.items():
        for variant in ["tradeseq"] + [f"de_peak_{s}_logp" for s in SCHEMES] + ["de_peak_fitted"]:
            t = sc_tab[variant].sort_values(["score", "waldStat_1"], ascending=False, kind="mergesort")
            weight = np.arange(len(t), 0, -1, dtype=float)
            for bench, gt in gts.items():
                for tf in g04b.CONTEXT_TFS[context]:
                    if tf not in set(gt["TF"]):
                        continue
                    e = pd.DataFrame({"TF": tf, "target": t["gene"].to_numpy(), "EdgeWeight": weight})
                    row, _ = b.evaluate_one(e[e["target"] != tf], gt, {tf}, valid, curve_points=50)
                    if row["k_pos"]:
                        pertf.append({"benchmark": bench, "context": context, "TF": tf, "variant": variant, "k_pos": row["k_pos"],
                                      "EPR": min(row["EPR"], pub_v2.EPR_CAP), "AUPRC_ratio": row["AUPRC_ratio"]})
    pertf = pd.DataFrame(pertf)
    pertf.to_csv(TABLES / "per_tf_by_peak_scheme.csv", index=False)
    five_b = {"K562": ["RELA", "NFKB1", "JUN", "CREB1", "ATF1", "ATF2"], "THP1": ["RELA", "NFKB1", "JUN", "ELK1", "CREB1", "ATF1", "ATF2", "REL", "IRF3"]}
    for (bench, context), g in pertf.groupby(["benchmark", "context"], sort=False):
        w = g.pivot(index="TF", columns="variant", values=["EPR", "AUPRC_ratio"]).round(2)
        print(f"\n== per-TF {bench} | {context} ==\n" + w.reindex([t for t in g04b.CONTEXT_TFS[context] if t in w.index]).to_string())
        print("median (5b set): " + g[g.TF.isin(five_b[context])].groupby("variant")[["EPR", "AUPRC_ratio"]].median().round(2).to_string().replace("\n", " | "))

    # -------------------------------------------------------------- global TF-balanced recall (04 engine)
    coexp = {}
    gt_tfs = set().union(*[set(g["TF"]) for g in gts.values()])
    for m in g04.COEXP:
        e = compare.c.load_beeline_method(m, valid)
        coexp[m] = e[e["TF"].isin(gt_tfs)]
    rows, curve_rows = [], []
    ks = g04.K_GRID.astype(float)
    for bench in g04.BENCHMARKS:
        gt = gts[bench]
        priors = {ctx: g04.chronopack_components(ctx, valid)[0] for ctx in g04b.CONTEXT}
        avail = [set(e["TF"]) for e in coexp.values()] + [set(p.index) for p in priors.values()]
        universe_tfs = sorted(set(gt["TF"]).intersection(*avail))
        targets_u = sorted(valid)
        n_universe = len(universe_tfs) * len(targets_u)
        gt_u = gt[gt["TF"].isin(universe_tfs) & gt["target"].isin(valid) & (gt["TF"] != gt["target"])].drop_duplicates(["TF", "target"])
        positives = set(zip(gt_u["TF"], gt_u["target"]))
        pos_by_tf = gt_u["TF"].value_counts().to_dict()
        for context, sc_tab in target_scores.items():
            signaling = [t for t in g04.CONTEXT_TFS[context] if t in universe_tfs and pos_by_tf.get(t, 0) >= g04.MIN_TARGETS]
            sets = {"min3": signaling, "min3_noRELA_NFKB1": [t for t in signaling if t not in g04.EXCLUDED_DOMINANT]}
            for variant in GLOBAL_VARIANTS:
                e = g04.outer_edges(priors[context], sc_tab[variant], universe_tfs, use_prior=True, use_target=True)
                e = e[e["target"].isin(valid)]
                rec = g04.Ranking(e, n_universe).recall(positives, {t: pos_by_tf.get(t, 0) for t in signaling}, signaling, ks)
                for set_name, tfl in sets.items():
                    curve = rec.loc[tfl].mean(axis=0)
                    rows.append(g04.summarise(bench, context, "all", set_name, variant, tfl, rec, curve, n_universe, rng))
                    for kk, v in curve.items():
                        curve_rows.append({"benchmark": bench, "context": context, "rewarded": set_name, "variant": variant, "K": kk, "R_bal": v})
    glob = pd.DataFrame(rows)
    glob.to_csv(TABLES / "global_calibration_summary.csv", index=False)
    curves = pd.DataFrame(curve_rows)
    curves.to_csv(TABLES / "global_calibration_curves.csv", index=False)
    cols = ["method", "logAUC", "logAUC_ci_lo", "logAUC_ci_hi", "logAUC_vs_random", "R@1000", "R@10000", f"n_TFs_recall>={g04.MEANINGFUL_RECALL:g}@10000"]
    for (bench, context, rew), g in glob.groupby(["benchmark", "context", "rewarded"], sort=False):
        print(f"\n== global {bench} | {context} | {rew} ==\n" + g[cols].round(4).to_string(index=False))

    make_figures(target_scores, peaks_all, curves, pertf)
    print(f"\nWrote {TABLES} and {FIGS}")


def make_figures(target_scores, peaks_all, curves, pertf) -> None:
    plt.rcParams.update({"pdf.fonttype": 42, "font.size": 8})
    # 1. Score calibration: score vs rank.
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    colors = {"tradeseq": "#111111", "de_minp": "#7570b3", "de_peak_pcdf": "#e6ab02", "de_peak_fitted": "#d95f02", "de_peak_observed": "#1b9e77"}
    labels = {"tradeseq": "tradeSeq 1-lfdr", "de_minp": "DE min-p 1-lfdr (04)", "de_peak_pcdf": "peak DE 1-lfdr (pc_df peak)",
              "de_peak_fitted": "peak DE 1-lfdr (fitted peak)", "de_peak_observed": "peak DE 1-lfdr (observed peak)"}
    for ax, (context, sc_tab) in zip(axes, target_scores.items()):
        for k, c in colors.items():
            s = np.sort(sc_tab[k]["score"].to_numpy())[::-1]
            ax.plot(np.arange(1, len(s) + 1), s, color=c, lw=1.4, label=labels[k])
        ax.set_xscale("log")
        ax.set_xlabel("Gene rank")
        ax.set_ylabel("Target score")
        ax.set_title(f"{context}: target-score calibration", loc="left")
        ax.grid(True, color="#eeeeee")
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGS / f"1_score_calibration.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 2. Global TF-balanced recall for calibration variants (TRRUST).
    vcolors = {"tradeseq": "#111111", "de_minp": "#7570b3", "de_peak": "#d95f02", "de_minp_qm": "#7570b3", "de_peak_qm": "#d95f02", "tradeseq_qm_de": "#111111"}
    vstyle = {"tradeseq": "-", "de_minp": "-", "de_peak": "-", "de_minp_qm": "--", "de_peak_qm": "--", "tradeseq_qm_de": ":"}
    fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharex=True)
    for i, bench in enumerate(g04.BENCHMARKS):
        for j, (context, rew) in enumerate([(c, r) for c in g04b.CONTEXT for r in ("min3", "min3_noRELA_NFKB1")]):
            ax = axes[i, j]
            d = curves[(curves.benchmark == bench) & (curves.context == context) & (curves.rewarded == rew)]
            for v, lab in GLOBAL_VARIANTS.items():
                dv = d[d.variant == v].sort_values("K")
                ax.plot(dv.K, dv.R_bal, color=vcolors[v], ls=vstyle[v], lw=1.6, label=lab)
            ax.set_xscale("log")
            ax.set_title(f"{bench} | {context} | {'all signaling TFs' if rew == 'min3' else 'without RELA/NFKB1'}", loc="left", fontsize=8)
            ax.grid(True, color="#eeeeee")
            if j == 0:
                ax.set_ylabel("TF-balanced recall")
            if i == 1:
                ax.set_xlabel("Global edge budget K")
            if i == 0 and j == 0:
                ax.legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGS / f"2_global_calibration_variants.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 3. Peak-time diagnostics: pc_df (converted) vs fitted peak minute.
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, (context, peaks) in zip(axes, peaks_all.items()):
        cell_type, mpt, conv, t_max = g04b.CONTEXT[context]
        pc = pd.read_csv(g04b.PC_DF)
        pc = pc[pc["cell_type"] == cell_type].set_index("gene")
        gs = [g for g in pc.index if g in peaks.index]
        x = pc.loc[gs, "t_peak"].map(conv)
        y = peaks.loc[gs, "fitted_peak_min"]
        capped = (pc.loc[gs, "t_peak"] >= 97 - 1e-9) if context == "K562" else pd.Series(False, index=gs)
        ax.scatter(x[~capped], y[~capped], s=4, alpha=0.35, color="#444444", label="pc_df genes")
        if capped.any():
            ax.scatter(x[capped], y[capped], s=6, alpha=0.6, color="#d95f02", label="K562 t_peak at the 97-min cap")
        ax.plot([0, t_max], [0, t_max], color="#c0392b", lw=0.8)
        ax.set_xlabel("pc_df t_peak converted to context minutes")
        ax.set_ylabel("Fitted-curve peak (screen fit, native axis), min")
        ax.set_title(f"{context}: peak-time conversion check", loc="left")
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGS / f"3_peak_time_check.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
