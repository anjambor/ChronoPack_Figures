#!/usr/bin/env python
"""Compare Fig4 v2 with v1. Writes output/comparison_to_v1.txt.

1. Rankings, panels A-B: the v2 manifest and positive-control ranks must equal v1's (fails if not).
2. Fit check: the quick fit reproduces the stored r9 final fit for the ligands in it.
3. Trajectories: per curve, the fitted mean (at the 12 sampled time points) vs the v1 curve
   (per-time-point mean of the log1p layer): Pearson r, and the peak minute of each.
4. Per method: top-10 ligands whose curve peaks inside the 0-45 min causal window, v2 vs v1.
"""
from __future__ import annotations

import anndata as ad
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
    OUTPUT_DIR,
    RECEIVER,
    RESPONSE_GENE,
    SENDER,
    TIME_KEY,
    V1_DIR,
    time_to_minute,
    to_dense,
)

V1_PREFIX = "fig4_composite_temporal_prioritization"
V2_PREFIX = "fig4_v2_composite_temporal_prioritization"
lines: list[str] = []


def say(text: str = "") -> None:
    print(text)
    lines.append(text)


def observed_means(adata, genes: list[str], cell_type: str) -> pd.DataFrame:
    mask = adata.obs[CELL_TYPE_KEY].astype(str).to_numpy() == cell_type
    sub = adata[mask, genes]
    x = np.asarray(to_dense(sub.layers[EXPR_LAYER]), dtype=float)
    df = pd.DataFrame(x, columns=genes)
    df["minute"] = time_to_minute(sub.obs[TIME_KEY]).to_numpy()
    return df.groupby("minute").mean().sort_index()


def main() -> None:
    failures = []
    say("== 1. Rankings and panels A-B: v2 vs v1 ==")
    for suffix in ["_manifest.csv", "_positive_control_ranks.csv"]:
        v2 = pd.read_csv(FIGURE_DIR / f"{V2_PREFIX}{suffix}")
        v1 = pd.read_csv(V1_DIR / f"{V1_PREFIX}{suffix}")
        if suffix == "_positive_control_ranks.csv":  # v2 adds tie-block columns (2026-11-16); compare v1's columns
            say("  positive-control tie blocks: " + "; ".join(
                f"{r.method} {r.rank_min:g}-{r.rank_max:g} (n={r.tie_n}, mid {r.mid_rank:g})" for r in v2.itertuples()))
            v2 = v2[v1.columns]
        try:
            pd.testing.assert_frame_equal(v2, v1, check_exact=False, rtol=1e-12, atol=1e-12)
            say(f"  {suffix[1:]}: identical ({len(v2)} rows)")
        except AssertionError as exc:
            failures.append(suffix)
            say(f"  {suffix[1:]}: DIFFERS\n{exc}")

    summary = pd.read_csv(FIT_SUMMARY_CSV)
    checked = summary.dropna(subset=["max_abs_diff_vs_r9_final_fit"])
    say(f"\n== 2. Quick fit vs stored r9 THP1 final fit (sce.1trajectory.nknots8.rds) ==")
    say(f"  {len(checked)} ligands in the stored fit; max |diff mu_log1p| = {checked['max_abs_diff_vs_r9_final_fit'].max():.3g}")
    nt = summary.loc[~summary["in_r9_tradeseq_inputs"].astype(bool)]
    say(f"  Not in the r9 THP1 tradeSeq inputs (fitted from unfiltered counts): "
        + ", ".join(f"{r.gene} ({int(r.thp1_total_counts)} counts)" for r in nt.itertuples()))

    curves = pd.read_csv(CURVES_CSV)
    manifest = pd.read_csv(FIGURE_DIR / f"{V2_PREFIX}_manifest.csv")
    adata = ad.read_h5ad(ADATA, backed="r")
    ligands = sorted(curves.loc[curves["cell_type"] == SENDER, "gene"].unique())
    obs = {SENDER: observed_means(adata, ligands, SENDER), RECEIVER: observed_means(adata, [RESPONSE_GENE], RECEIVER)}

    say("\n== 3. Trajectories: fitted (v2) vs per-time-point mean log1p (v1) ==")
    say(f"  {'gene':8s} {'cell':5s} {'r (12 points)':>13s} {'peak v2':>8s} {'peak v1':>8s}  r9 tested")
    rows = []
    for (ct, gene), g in curves.groupby(["cell_type", "gene"]):
        fit = g.set_index("t")["mu_log1p"]
        o = obs[ct][gene]
        f12 = fit.reindex(o.index.astype(float))
        r = float(np.corrcoef(f12, o)[0, 1]) if o.std() > 0 and f12.std() > 0 else np.nan
        tested = summary.set_index("gene")["in_r9_tradeseq_inputs"].get(gene, True)
        rows.append({"cell_type": ct, "gene": gene, "r": r, "peak_v2": float(fit.idxmax()), "peak_v1": float(o.idxmax())})
        say(f"  {gene:8s} {ct:5s} {r:13.3f} {fit.idxmax():8.0f} {o.idxmax():8.0f}  {'yes' if tested else 'no'}")
    traj = pd.DataFrame(rows)
    traj.to_csv(OUTPUT_DIR / "chronopack" / "fitted_vs_observed_curves.csv", index=False)
    say(f"  median r = {traj['r'].median():.3f}; r < 0.8: {', '.join(traj.loc[traj['r'] < 0.8, 'gene']) or 'none'}")

    say(f"\n== 4. Top-10 ligands peaking in the causal window ({CAUSAL_WINDOW[0]:.0f}-{CAUSAL_WINDOW[1]:.0f} min), v2 (fit) vs v1 (means) ==")
    peaks = traj.loc[traj["cell_type"] == SENDER].set_index("gene")
    for spec in CURVE_METHODS:
        top = manifest.loc[manifest["plot_method"] == spec["label"], "ligand"].astype(str)
        in_win = lambda col: sorted(g for g in top if CAUSAL_WINDOW[0] <= peaks.loc[g, col] <= CAUSAL_WINDOW[1])
        v2, v1 = in_win("peak_v2"), in_win("peak_v1")
        say(f"  {spec['label']:10s} v2 {len(v2):2d} ({', '.join(v2) or '-'})   v1 {len(v1):2d} ({', '.join(v1) or '-'})")

    (OUTPUT_DIR / "comparison_to_v1.txt").write_text("\n".join(lines) + "\n")
    if failures:
        raise SystemExit(f"v2 rankings differ from v1: {failures}")
    print("COMPARE_OK")


if __name__ == "__main__":
    main()
