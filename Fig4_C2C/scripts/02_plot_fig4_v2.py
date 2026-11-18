#!/usr/bin/env python
"""Fig4 v2: composite temporal prioritization of THP1 -> K562 ligands (top 10 per method).

Fig4_C2C/scripts/07_plot_fig4.py with one change: the trajectory panels (right) draw the
ChronoPack NB-GAM fits of step 01 (THP1 ligands, K562 NFKBIA) instead of per-time-point means
of the r9 adata log1p layer. Rankings come from v1's tables (input/v1/); panels A and B and the
manifest (early log2FC etc. from the log1p layer) are computed exactly as in v1.
"""
from __future__ import annotations

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig4_v2_common import (
    ADATA,
    CAUSAL_WINDOW,
    CELL_TYPE_KEY,
    CURVE_METHODS,
    CURVES_CSV,
    EXPR_LAYER,
    FIGURE_DIR,
    FIT_SUMMARY_CSV,
    K562_FIT_DIR,
    METHODS,
    NKNOTS,
    POSITIVE_CONTROL_PAIRS,
    RECEIVER,
    RESPONSE_GENE,
    SENDER,
    THP1_FIT_DIR,
    TIME_KEY,
    TRASIG_N_PERMS,
    build_top_manifest,
    method_rows,
    time_to_minute,
    to_dense,
)

OUT_PREFIX = FIGURE_DIR / "fig4_v2_composite_temporal_prioritization"
RESPONSE_GENES = [RESPONSE_GENE]
LOG2FC_PSEUDOCOUNT = 0.01
CURVE_TOP_N = 10
EARLY_LOG2FC_BASE_MIN = 0
EARLY_LOG2FC_COMPARE_MIN = 45


def mean_signal_by_time(adata, genes: list[str], obs_mask=None) -> pd.Series:
    genes = [gene for gene in genes if gene in adata.var_names]
    if not genes:
        return pd.Series(dtype=float)
    sub = adata[:, genes] if obs_mask is None else adata[obs_mask, genes]
    values = np.asarray(to_dense(sub.layers[EXPR_LAYER]), dtype=float).mean(axis=1)
    df = pd.DataFrame({"minute": time_to_minute(sub.obs[TIME_KEY]).to_numpy(), "expr": values})
    return df.dropna().groupby("minute")["expr"].mean().sort_index()


def minmax(series: pd.Series) -> pd.Series:
    if series.empty:
        return series
    values = series.astype(float)
    lo, hi = float(values.min()), float(values.max())
    if not np.isfinite(hi - lo) or hi == lo:
        return values * 0.0
    return (values - lo) / (hi - lo)


def response_onset_min(response_norm: pd.Series, threshold: float = 0.5) -> float:
    post = response_norm.loc[response_norm.index > 0]
    hit = post.loc[post >= threshold]
    return float(hit.index.min()) if not hit.empty else float(response_norm.index.max())


def expression_at_minute(series: pd.Series, minute: float) -> float:
    if series.empty:
        return np.nan
    matches = series.loc[np.isclose(series.index.astype(float), minute)]
    return float(matches.iloc[0]) if not matches.empty else np.nan


def ligand_metrics(manifest: pd.DataFrame, ligand_series: dict[str, pd.Series], response_norm: pd.Series) -> pd.DataFrame:
    onset = response_onset_min(response_norm)
    rows = []
    for _, row in manifest.iterrows():
        raw_series = ligand_series.get(str(row["ligand"]), pd.Series(dtype=float))
        series = minmax(raw_series)
        expr_base = expression_at_minute(raw_series, EARLY_LOG2FC_BASE_MIN)
        expr_compare = expression_at_minute(raw_series, EARLY_LOG2FC_COMPARE_MIN)
        early_log2fc = (
            float(np.log2((expr_compare + LOG2FC_PSEUDOCOUNT) / (expr_base + LOG2FC_PSEUDOCOUNT)))
            if np.isfinite(expr_base) and np.isfinite(expr_compare)
            else np.nan
        )
        if series.empty:
            peak_min = amp = baseline_mean = causal_mean = causal_enrichment = np.nan
        else:
            peak_min = float(series.idxmax())
            amp = float(series.max() - series.min())
            baseline_mean = float(series.loc[(series.index >= 0) & (series.index <= 15)].mean())
            causal_mean = float(series.loc[(series.index >= CAUSAL_WINDOW[0]) & (series.index <= CAUSAL_WINDOW[1])].mean())
            causal_enrichment = causal_mean - baseline_mean
        plausible = bool(np.isfinite(peak_min) and CAUSAL_WINDOW[0] <= peak_min <= min(CAUSAL_WINDOW[1], onset))
        rows.append(
            {
                **row.to_dict(),
                "ligand_peak_min": peak_min,
                "dynamic_amplitude": amp,
                "baseline_mean_norm": baseline_mean,
                "causal_window_mean_norm": causal_mean,
                "causal_window_enrichment": causal_enrichment,
                "expression_t0": expression_at_minute(raw_series, 0),
                "expression_t15": expression_at_minute(raw_series, 15),
                "expression_t45": expression_at_minute(raw_series, 45),
                "early_log2fc_45_vs_0": early_log2fc,
                "plausible_window_peak": plausible,
            }
        )
    return pd.DataFrame(rows)


def positive_control_ranks() -> pd.DataFrame:
    """Rank of the positive control per method. `rank` is the method's own (v1) rank; rows tied with it on the
    method's score form a tie block (rank_min-rank_max, tie_n rows), and mid_rank is the block's average rank
    (panel A bar length; user request 2026-11-16: CellChat's TNF pair is one of 28 pairs tied at lr_probs = 0)."""
    rows = []
    for spec in METHODS:
        data = method_rows(spec["label"]).sort_values("rank", kind="stable")
        if spec["label"] == "NicheNet":
            hit = data.loc[data["ligand"].astype(str) == "TNF"].head(1)
            target = "TNF ligand"
        else:
            hit = data.loc[data["interaction_name"].astype(str).isin(POSITIVE_CONTROL_PAIRS)]
            hit = hit.loc[hit["rank"] == hit["rank"].min()]
            target = "best TNF/TNFR pair"
        rank = rank_min = rank_max = mid_rank = np.nan
        tie_n = 0
        if not hit.empty:
            rank = float(hit["rank"].iloc[0])
            score = float(hit["score"].iloc[0])
            block = data.loc[np.isclose(data["score"].astype(float), score, rtol=1e-12, atol=1e-15), "rank"].astype(float)
            rank_min, rank_max, tie_n = float(block.min()), float(block.max()), int(len(block))
            assert rank_min <= rank <= rank_max and rank_max - rank_min + 1 == tie_n, (spec["label"], block.tolist())
            mid_rank = (rank_min + rank_max) / 2
        rows.append({"method": spec["label"], "rank": rank, "target": target, "rank_min": rank_min, "rank_max": rank_max,
                     "tie_n": tie_n, "mid_rank": mid_rank, "color": spec["color"]})
    return pd.DataFrame(rows)


def plot_rank_panel(ax, rank_df: pd.DataFrame) -> None:
    y = np.arange(len(rank_df))
    ranks = rank_df["mid_rank"].astype(float).fillna(100)  # tied blocks: bar at the block's mid-rank
    ax.barh(y, ranks, color=rank_df["color"].tolist(), alpha=0.9)
    for yi, row in zip(y, rank_df.itertuples(index=False)):
        if pd.isna(row.mid_rank):
            text = "not in ranked output"
        elif row.tie_n > 1:
            text = f"rank {int(row.rank_min)}-{int(row.rank_max)} (tie)"
        else:
            text = f"rank {int(row.rank)}"
        ax.text((0 if pd.isna(row.mid_rank) else row.mid_rank) + 1.0, yi, text, va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(rank_df["method"], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Rank (lower is better)", fontsize=9)
    ax.set_title("A. Positive-control recovery", fontsize=11, loc="left")
    ax.set_xlim(0, max(25, float(ranks.max()) + 16))  # room for the "(tie)" labels
    ax.grid(True, axis="x", color="#e5e5e5", lw=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_plausibility_panel(ax, metrics: pd.DataFrame) -> None:
    labels = [spec["label"] for spec in METHODS]
    colors = [spec["color"] for spec in METHODS]
    values = [metrics.loc[metrics["plot_method"] == label, "early_log2fc_45_vs_0"].dropna().to_numpy(dtype=float) for label in labels]
    positions = np.arange(len(labels), 0, -1)
    violin = ax.violinplot(values, positions=positions, vert=False, widths=0.72, showmedians=True, showextrema=False)
    for body, color in zip(violin["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.28)
        body.set_linewidth(1.0)
    violin["cmedians"].set_color("0.2")
    violin["cmedians"].set_linewidth(1.2)

    rng = np.random.default_rng(11)
    for pos, color, vals in zip(positions, colors, values):
        jitter = rng.normal(0, 0.035, size=len(vals))
        ax.scatter(vals, np.full(len(vals), pos) + jitter, s=24, color=color, edgecolor="white", linewidth=0.5, zorder=3)

    all_values = np.concatenate([v for v in values if len(v)])
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlim(min(-0.4, float(np.nanmin(all_values)) - 0.2), max(4.8, float(np.nanmax(all_values)) + 0.35))
    ax.set_xlabel("Early log2FC (45 min vs. 0 min)", fontsize=9)
    ax.set_title("B. Early ligand induction", fontsize=11, loc="left")
    ax.grid(True, axis="x", color="#e5e5e5", lw=0.8)
    ax.axvline(0, color="#e5e5e5", lw=0.8, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def fitted_series(curves: pd.DataFrame, cell_type: str) -> dict[str, pd.Series]:
    sub = curves.loc[curves["cell_type"] == cell_type]
    return {gene: g.set_index("t")["mu_log1p"].sort_index() for gene, g in sub.groupby("gene")}


def plot_curve_panel(ax, rows: pd.DataFrame, ligand_series: dict[str, pd.Series], response_expr: pd.Series, spec: dict, not_tested: set[str]) -> None:
    cmap = plt.get_cmap("tab10")
    for _, row in rows.sort_values("ligand_rank").head(CURVE_TOP_N).iterrows():
        ligand = str(row["ligand"])
        series = ligand_series.get(ligand, pd.Series(dtype=float))
        if series.empty:
            raise ValueError(f"No fitted curve for {ligand}; rerun step 01 (REFIT=1).")
        rank = int(row["ligand_rank"])
        mark = "\u2020" if ligand in not_tested else ""
        ax.plot(series.index, series, lw=1.35, color=cmap((rank - 1) % 10), alpha=0.88, label=f"{rank}. {ligand}{mark}")
    ax.plot(response_expr.index, response_expr, color="black", lw=2.8, label=f"{RECEIVER} NFKBIA")
    ax.axvspan(CAUSAL_WINDOW[0], CAUSAL_WINDOW[1], color=spec["color"], alpha=0.08, lw=0)
    ax.set_xlim(0, 165)
    ax.set_ylim(bottom=-0.02)
    ax.set_xticks([0, 30, 60, 90, 120, 150])
    ax.set_title(spec["label"], fontsize=10.5, loc="left")
    ax.grid(True, color="#e7e7e7", lw=0.8)
    ax.tick_params(labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5), borderaxespad=0)


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    manifest = build_top_manifest()
    adata = ad.read_h5ad(ADATA, backed="r")
    cell_types = adata.obs[CELL_TYPE_KEY].astype(str).to_numpy()
    response_expr = mean_signal_by_time(adata, RESPONSE_GENES, obs_mask=cell_types == RECEIVER)
    response_norm = minmax(response_expr)
    ligands = sorted(set(manifest["ligand"].astype(str)))
    ligand_series = {ligand: mean_signal_by_time(adata, [ligand], obs_mask=cell_types == SENDER) for ligand in ligands}
    metrics = ligand_metrics(manifest, ligand_series, response_norm)
    metrics.to_csv(f"{OUT_PREFIX}_manifest.csv", index=False)

    curves = pd.read_csv(CURVES_CSV)
    fit_ligands = fitted_series(curves, SENDER)
    fit_response = fitted_series(curves, RECEIVER)[RESPONSE_GENE]
    summary = pd.read_csv(FIT_SUMMARY_CSV)
    not_tested = set(summary.loc[~summary["in_r9_tradeseq_inputs"].astype(bool), "gene"].astype(str))
    ref_lib = {ct: float(np.exp(np.log(pd.read_csv(d / "cell_libsize_full.csv")["libsize_full"].astype(float)).mean()))
               for ct, d in [(SENDER, THP1_FIT_DIR), (RECEIVER, K562_FIT_DIR)]}
    rank_df = positive_control_ranks()
    rank_df.drop(columns=["color"]).to_csv(f"{OUT_PREFIX}_positive_control_ranks.csv", index=False)

    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(13.4, 9.6), dpi=220)
    outer = fig.add_gridspec(1, 2, width_ratios=[1.02, 1.11], wspace=0.25)
    left = outer[0].subgridspec(2, 1, hspace=0.55)
    right = outer[1].subgridspec(len(CURVE_METHODS), 1, hspace=0.5)

    plot_rank_panel(fig.add_subplot(left[0]), rank_df)
    plot_plausibility_panel(fig.add_subplot(left[1]), metrics)
    for i, spec in enumerate(CURVE_METHODS):
        ax = fig.add_subplot(right[i])
        plot_curve_panel(ax, metrics.loc[metrics["plot_method"] == spec["label"]].copy(), fit_ligands, fit_response, spec, not_tested)
        ax.set_xlabel("Time (min)" if i == len(CURVE_METHODS) - 1 else "", fontsize=9)
        if i == len(CURVE_METHODS) // 2:
            ax.set_ylabel("Fitted expression (log1p)", fontsize=9)

    fig.suptitle("ChronoPack prioritizes dynamically timed ligand-receptor signals before K562 NF-kB activation", fontsize=13.5, y=0.982)
    fig.text(
        0.5,
        0.015,
        "ChronoPack r9 data. NicheNet uses THP1 ligand detection over 0-75 min and K562b peak-response targets (peak log2FC > 0.5, p < 0.05) "
        f"ranked by corrected AUPR; TraSig uses {TRASIG_N_PERMS:,} permutations. Part B shows top-10 ligand log2FC distributions from 0 to 45 min "
        "with pseudocount 0.01.\n"
        f"Trajectories show ChronoPack NB-GAM fits ({NKNOTS} knots; log1p fitted mean count at the reference library size, "
        f"{SENDER} {ref_lib[SENDER]:,.0f} / {RECEIVER} {ref_lib[RECEIVER]:,.0f} UMIs) for the top 10 ligands ({SENDER}) and {RECEIVER} {RESPONSE_GENE}, "
        "excluding CellChat. \u2020 not in the r9 THP1 tradeSeq inputs (low-count or technical filter); fitted from unfiltered counts.",
        ha="center",
        fontsize=8.5,
        color="0.25",
    )
    fig.subplots_adjust(left=0.075, right=0.905, top=0.925, bottom=0.125, wspace=0.25)
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{OUT_PREFIX}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(rank_df.drop(columns=["color"]).to_string(index=False))
    for ext in ["png", "pdf", "svg"]:
        print(f"{OUT_PREFIX}.{ext}")


if __name__ == "__main__":
    main()
