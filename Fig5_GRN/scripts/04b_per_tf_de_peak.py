#!/usr/bin/env python
"""Per-TF GRN performance with ChronoPack target scores from peak-time differential expression
(side analysis; Fig5 outputs are not changed).

ChronoPack's edge score is TF prior x target score. The TF prior is constant within a TF, so per-TF
EPR / AUPRCr depend only on how targets are ordered; this script swaps the target score:

  tradeseq_lfdr   1 - lfdr_bum of the tradeSeq association p-value (current Fig5 ChronoPack; Wald tie order)
  peak_de_logp    -log10 p of a Wilcoxon test, cells at the gene's peak time vs baseline (|z| tie order)
  peak_de_lfdr    1 - lfdr_bum of that p (ChronoPack's own BUM calibration, cp.tl.grn.score_targets)
  peak_abs_lfc    |log2 fold change| at the peak time vs baseline (pseudobulk means; |z| tie order)

Peak time per gene: for genes in the r9 gene-trajectory table (pc_df.new.csv, identical to
pc_df.updated.csv for K562a and THP1; tradeSeq p < 0.05 genes), its t_peak mapped back from the aligned
clustering axis to the context's own minutes (K562a: minutes as is, capped at 97 by the shared window;
THP1: t_peak + 10 / 0.65) and then to the nearest sampled time point. For the other tested genes, the
time point with the largest |log2FC| of mean normalized expression vs baseline.

Cells, genes and counts are those of the r9 tradeSeq fit (input/chronopack/<context>/tradeseq_fit),
normalized to 10,000 per cell (full library size) and log1p-transformed (as in 04a). Evaluation: the
Fig5 per-TF evaluation (b.evaluate_one, valid targets, publication ground truths) for the Fig5
signaling TFs. Writes output/per_tf_de_peak/.
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
import scanpy as sc
import anndata as ad

import chronopack as cp
import benchmark_grn_systems as b
import generate_system_broad_filtered_publication_figures_v2 as pub_v2
from plot_regulator_prioritization_v2 import K562_EVAL_TFS, THP1_EVAL_TFS

OUT = fig5_paths.OUTPUT_DIR / "per_tf_de_peak"
PROTEIN_CODING = fig5_paths.CHRONOPACK_INPUT_DIR / "protein_coding_genes.csv"
PC_DF = Path("/new-stg/home/alex/tmp/ChronoPack/verification_runs/full_pipeline_analyses_20260930_r9/data/gene_clustering/pc_df.new.csv")
CONTEXT = {  # context: (pc_df cell_type, minutes per time tag, t_peak -> context minutes, max context minute)
    "K562": ("K562a", 10.0, lambda t: t, 110.0),
    "THP1": ("THP1", 15.0, lambda t: t + 10.0 / 0.65, 165.0),
}
CONTEXT_TFS = {"K562": list(K562_EVAL_TFS), "THP1": list(THP1_EVAL_TFS)}
LFC_PSEUDOCOUNT = 0.1  # on the CP10k scale, for choosing the observed peak of genes outside pc_df
VARIANTS = ["tradeseq_lfdr", "peak_de_logp", "peak_de_lfdr", "peak_abs_lfc"]
LABELS = {"tradeseq_lfdr": "tradeSeq 1-lfdr (Fig5)", "peak_de_logp": "Peak DE -log10 p", "peak_de_lfdr": "Peak DE 1-lfdr", "peak_abs_lfc": "Peak |log2FC|"}


def load_fit(context: str) -> ad.AnnData:
    fit = fig5_paths.CHRONOPACK_INPUT_DIR / context / "tradeseq_fit"
    counts = pd.read_csv(fit / "counts_genes_x_cells.csv.gz", index_col=0)
    times = pd.read_csv(fit / "cell_times.csv.gz", index_col="cell_id")["time"]
    lib = pd.read_csv(fit / "cell_libsize_full.csv", index_col="cell_id")["libsize_full"]
    cells = counts.columns.astype(str)
    cp10k = counts.T.to_numpy(dtype=np.float32) / lib.reindex(cells).to_numpy(float)[:, None] * 1e4
    a = ad.AnnData(X=np.log1p(cp10k), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=counts.index.astype(str)))
    a.layers["cp10k"] = cp10k
    tags = sorted(times.unique())
    a.obs["time"] = pd.Categorical(times.reindex(cells).astype(int).astype(str), categories=[str(t) for t in tags])
    return a


def peak_tags(context: str, a: ad.AnnData) -> pd.DataFrame:
    cell_type, minutes_per_tag, to_minutes, t_max = CONTEXT[context]
    tags = [int(t) for t in a.obs["time"].cat.categories]
    means = pd.DataFrame({t: np.asarray(a.layers["cp10k"][(a.obs["time"] == str(t)).to_numpy()].mean(axis=0)).ravel() for t in tags}, index=a.var_names)
    lfc = np.log2((means + LFC_PSEUDOCOUNT).div(means[tags[0]] + LFC_PSEUDOCOUNT, axis=0))
    later = lfc[tags[1:]]
    observed_tag = later.abs().idxmax(axis=1).astype(int)
    pc = pd.read_csv(PC_DF)
    pc = pc[pc["cell_type"] == cell_type].set_index("gene")
    minutes = pc["t_peak"].map(to_minutes).clip(lower=0, upper=t_max)
    pc_tag = (np.round(minutes / minutes_per_tag).astype(int) + tags[0]).clip(lower=tags[1], upper=tags[-1])
    out = pd.DataFrame({"observed_peak_tag": observed_tag})
    out["pc_df_peak_tag"] = pc_tag.reindex(out.index)
    out["pc_df_peak_log2fc"] = pc["peak_log2fc"].reindex(out.index)
    out["peak_tag"] = out["pc_df_peak_tag"].fillna(out["observed_peak_tag"]).astype(int)
    out["peak_source"] = np.where(out["pc_df_peak_tag"].notna(), "pc_df", "observed")
    out["log2fc_at_peak"] = [lfc.loc[g, t] for g, t in zip(out.index, out["peak_tag"])]
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    valid = b.load_valid_targets()
    gts = pub_v2.load_reference_ground_truths(valid)
    published = pd.read_csv(fig5_paths.PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv")

    rows = []
    for context in CONTEXT:
        a = load_fit(context)
        peaks = peak_tags(context, a)
        tags = list(a.obs["time"].cat.categories)
        sc.tl.rank_genes_groups(a, groupby="time", groups=tags[1:], reference=tags[0], method="wilcoxon", n_genes=a.n_vars)
        res = a.uns["rank_genes_groups"]
        p_tab = pd.DataFrame({int(t): pd.Series(res["pvals"][t], index=res["names"][t]) for t in tags[1:]}).reindex(a.var_names)
        z_tab = pd.DataFrame({int(t): pd.Series(res["scores"][t], index=res["names"][t]) for t in tags[1:]}).reindex(a.var_names)
        peaks["p_at_peak"] = [p_tab.loc[g, t] for g, t in zip(peaks.index, peaks["peak_tag"])]
        peaks["z_at_peak"] = [z_tab.loc[g, t] for g, t in zip(peaks.index, peaks["peak_tag"])]
        pcg = peaks["peak_source"] == "pc_df"
        agree = (peaks.loc[pcg, "pc_df_peak_tag"] == peaks.loc[pcg, "observed_peak_tag"]).mean()
        within1 = ((peaks.loc[pcg, "pc_df_peak_tag"] - peaks.loc[pcg, "observed_peak_tag"]).abs() <= 1).mean()
        sign = (np.sign(peaks.loc[pcg, "pc_df_peak_log2fc"]) == np.sign(peaks.loc[pcg, "log2fc_at_peak"])).mean()
        print(f"{context}: {a.n_obs} cells, {a.n_vars} genes; pc_df genes {int(pcg.sum())}: pc_df peak = observed peak {agree:.1%} "
              f"(within one time point {within1:.1%}); sign of fold change agrees {sign:.1%}; "
              f"peak-DE p < 0.05: {(peaks.p_at_peak < 0.05).sum()}", flush=True)
        peaks.to_csv(OUT / f"peak_de_{context}.csv")

        assoc = OUT / f"peak_de_association_{context}.csv"
        pd.DataFrame({"pvalue_1": peaks["p_at_peak"], "waldStat_1": peaks["z_at_peak"].abs()}, index=peaks.index.rename("gene")).to_csv(assoc)
        de_lfdr = cp.tl.grn.score_targets(assoc, valid_targets=PROTEIN_CODING).scores.set_index("gene")
        ts = pd.read_csv(fig5_paths.CHRONOPACK_OUTPUT_DIR / context / "target_scores.scores.csv").set_index("gene")

        genes = sorted(set(ts.index) & set(peaks.index) & valid)
        tie_z = peaks["z_at_peak"].abs().reindex(genes)
        scores = {
            "tradeseq_lfdr": (1 - ts["lfdr_bum"].reindex(genes), ts["waldStat_1"].reindex(genes)),
            "peak_de_logp": (-np.log10(peaks["p_at_peak"].reindex(genes).clip(lower=1e-300)), tie_z),
            "peak_de_lfdr": (1 - de_lfdr["lfdr_bum"].reindex(genes), tie_z),
            "peak_abs_lfc": (peaks["log2fc_at_peak"].abs().reindex(genes), tie_z),
        }
        for variant, (primary, secondary) in scores.items():
            order = pd.DataFrame({"target": genes, "s1": primary.to_numpy(), "s2": secondary.fillna(-np.inf).to_numpy()})
            order = order[np.isfinite(order["s1"]) & (order["s1"] > 0)].sort_values(["s1", "s2"], ascending=False, kind="mergesort")
            weight = np.arange(len(order), 0, -1, dtype=float)  # strictly decreasing: s1 then s2 decide the order
            for bench, gt in gts.items():
                for tf in CONTEXT_TFS[context]:
                    if tf not in set(gt["TF"]):
                        continue
                    edges = pd.DataFrame({"TF": tf, "target": order["target"].to_numpy(), "EdgeWeight": weight})
                    edges = edges[edges["target"] != tf]
                    row, _ = b.evaluate_one(edges, gt, {tf}, valid, curve_points=50)
                    if not row["k_pos"]:
                        continue
                    rows.append({"benchmark": bench, "context": context, "TF": tf, "variant": variant, "k_pos": row["k_pos"],
                                 "n_pred": row["n_pred"], "TP_at_k": row["TP_at_k"], "EPR": min(row["EPR"], pub_v2.EPR_CAP),
                                 "EPR_raw": row["EPR"], "AUPRC_ratio": row["AUPRC_ratio"]})
    res = pd.DataFrame(rows)

    # Check: the tradeseq_lfdr variant must reproduce the published per-TF ChronoPack metrics.
    pub = published[published["method"].isin(["ChronoPack_K562", "ChronoPack_THP1"])].copy()
    pub["context"] = pub["method"].str.replace("ChronoPack_", "", regex=False)
    chk = res[res.variant == "tradeseq_lfdr"].merge(pub[["benchmark", "context", "TF", "EPR", "AUPRC_ratio"]], on=["benchmark", "context", "TF"], suffixes=("", "_pub"))
    print(f"[check] tradeSeq variant vs published ({len(chk)} rows): max |dEPR| {np.abs(chk.EPR - chk.EPR_pub).max():.2e}, "
          f"max |dAUPRCr| {np.abs(chk.AUPRC_ratio - chk.AUPRC_ratio_pub).max():.2e}", flush=True)

    coexp = published[published["method"].isin(["GENIE3", "GRNBOOST2", "PPCOR", "SINCERITIES"])]
    best = coexp.groupby(["benchmark", "TF"]).agg(best_coexp_EPR=("EPR", "max"), best_coexp_AUPRCr=("AUPRC_ratio", "max")).reset_index()
    res = res.merge(best, on=["benchmark", "TF"], how="left")
    res.to_csv(OUT / "per_tf_metrics.csv", index=False)
    pd.set_option("display.width", 250)
    for (bench, context), g in res.groupby(["benchmark", "context"], sort=False):
        w = g.pivot(index="TF", columns="variant", values=["EPR", "AUPRC_ratio"]).reindex(columns=VARIANTS, level=1)
        extra = g.drop_duplicates("TF").set_index("TF")[["k_pos", "best_coexp_EPR", "best_coexp_AUPRCr"]]
        print(f"\n== {bench} | {context} ==")
        print(pd.concat([extra, w.round(2)], axis=1).reindex([t for t in CONTEXT_TFS[context] if t in extra.index]).to_string())
        med = g.groupby("variant")[["EPR", "AUPRC_ratio"]].median().reindex(VARIANTS).round(2)
        print("medians:\n" + med.to_string())
    plot(res)


def plot(res: pd.DataFrame) -> None:
    plt.rcParams.update({"pdf.fonttype": 42, "font.size": 8})
    combos = list(res.groupby(["benchmark", "context"], sort=False).groups)
    fig, axes = plt.subplots(len(combos), 2, figsize=(10, 2.6 * len(combos)))
    colors = {"tradeseq_lfdr": "#111111", "peak_de_logp": "#d95f02", "peak_de_lfdr": "#7570b3", "peak_abs_lfc": "#1b9e77"}
    for i, (bench, context) in enumerate(combos):
        g = res[(res.benchmark == bench) & (res.context == context)]
        tfs = [t for t in CONTEXT_TFS[context] if t in set(g.TF)]
        for j, metric in enumerate(["EPR", "AUPRC_ratio"]):
            ax = axes[i, j]
            x = np.arange(len(tfs))
            for k, v in enumerate(VARIANTS):
                vals = g[g.variant == v].set_index("TF")[metric].reindex(tfs)
                ax.bar(x + (k - 1.5) * 0.2, vals, width=0.2, color=colors[v], label=LABELS[v])
            best = g.drop_duplicates("TF").set_index("TF")[f"best_coexp_{'EPR' if metric == 'EPR' else 'AUPRCr'}"].reindex(tfs)
            ax.scatter(x, best, marker="_", s=260, color="#c0392b", zorder=5, label="best co-expression")
            ax.axhline(1, color="#999999", lw=0.8, ls="--")
            ax.set_xticks(x)
            ax.set_xticklabels(tfs)
            ax.set_title(f"{bench} | {context} | {'EPR (capped at 25)' if metric == 'EPR' else 'AUPRCr'}", loc="left")
            ax.spines[["top", "right"]].set_visible(False)
            if i == 0 and j == 1:
                ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"per_tf_de_peak.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
