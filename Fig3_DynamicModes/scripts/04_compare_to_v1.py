#!/usr/bin/env python
"""Compare the v2 (ChronoPack, r9) downsampling figure with v1 (old synthetic run) and sanity-check
the v2 reference fit. Writes output/comparison_to_v1.txt; reports differences, never fails on them (except check 5).

1. Per-mode n and mean +/- SE absolute offsets (bars), v2 vs v1.
2. Selected genes (q < 0.05 filter): shared / only-v1 / only-v2.
3. Per-gene onset/peak for the selected genes, v2 vs v1 (v1 times shifted by -1 to true minutes):
   Pearson r and mean absolute difference, for the reference and the downsampled fits.
4. Reference-fit check: v2 K562_trunc curves (0-90 min fit) vs r9's full 0-110 min K562a fit,
   restricted to 0-90 min, for the selected genes (r9 fit predicted with cp.tl.predict_gene_dynamics).
5. Broad/fine modes figure (step 00, added 2026-11-16): its CSVs must equal v1's and its PNG pixels are compared
   (v1 outputs and v1's comparison with the original Fig3 copied to input/v1/modes/). This check fails the step.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import chronopack as cp

FIG_DIR = Path(__file__).resolve().parents[1]
OUT = FIG_DIR / "output"
V1 = FIG_DIR / "input" / "v1"
STEM = "publication_timing_gene_composite_dynamic_modes_v10_original_dynamic_modes_v3_p005_fcgt05_all_assigned_combined_modes_with_transient_q005_no_ybreak_vector"
lines: list[str] = []


def say(text: str = "") -> None:
    print(text)
    lines.append(text)


def main() -> None:
    v2s, v1s = pd.read_csv(OUT / f"{STEM}_v2_mode_summary.csv"), pd.read_csv(V1 / f"{STEM}_mode_summary.csv")
    m = v2s.merge(v1s, on=["event", "dynamic_mode"], suffixes=("_v2", "_v1"))
    say("== Mean absolute offset per mode (min): v2 (ChronoPack, r9) vs v1 (old synthetic run) ==")
    for r in m.itertuples(index=False):
        say(f"  {r.event:5s} {r.dynamic_mode:15s} n={r.n_v2:4d} {r.mean_abs_offset_min_v2:6.2f} +/- {r.se_v2:4.2f}"
            f"   vs   n={r.n_v1:4d} {r.mean_abs_offset_min_v1:6.2f} +/- {r.se_v1:4.2f}")

    g2, g1 = pd.read_csv(OUT / f"{STEM}_v2_genes.csv"), pd.read_csv(V1 / f"{STEM}_genes.csv")
    s2, s1 = set(g2["gene"]), set(g1["gene"])
    say(f"\n== Selected genes: v2 {len(s2)}, v1 {len(s1)}, shared {len(s2 & s1)}, only v2 {len(s2 - s1)}, only v1 {len(s1 - s2)} ==")

    d2 = pd.read_csv(OUT / "chronopack" / "trajectory_descriptors.csv", index_col=0)
    d1 = pd.read_csv(V1 / "synthetic_ds_trajectory_descriptors.csv", index_col=0)
    d1 = d1.assign(t_onset=d1["t_onset"] - 1, t_peak=d1["t_peak"] - 1)  # v1 grid index -> true minutes
    genes = sorted(s2 & s1)
    say("\n== Per-gene timing of the selected genes, v2 vs v1 (true minutes) ==")
    for run in ["K562_trunc", "K562_ds"]:
        a = d2[d2["cell_type"].eq(run)].set_index("gene").reindex(genes)
        b = d1[d1["cell_type"].eq(run)].set_index("gene").reindex(genes)
        for col in ["t_onset", "t_peak"]:
            ok = a[col].notna() & b[col].notna()
            r = np.corrcoef(a.loc[ok, col], b.loc[ok, col])[0, 1]
            say(f"  {run:10s} {col:7s}: n={int(ok.sum())}, Pearson r={r:.3f}, mean |v2 - v1|={np.mean(np.abs(a.loc[ok, col] - b.loc[ok, col])):.2f} min")

    say("\n== Reference-fit check: v2 K562_trunc (0-90 min fit) vs r9 full 0-110 min fit, t <= 90 ==")
    curves = pd.read_csv(OUT / "chronopack" / "curves.csv")
    trunc = curves[curves["run"].eq("K562_trunc") & curves["gene"].isin(s2)]
    full = cp.tl.predict_gene_dynamics(FIG_DIR / "input" / "r9_full_fit", sorted(s2), n_points=111, time_max=110,
                                       sce_file="sce.1trajectory.nknots8.rds")
    full = full[full["t"].round(6) <= 90].assign(t=lambda d: d["t"].round(6))
    j = trunc.assign(t=trunc["t"].round(6)).merge(full, on=["gene", "t"], suffixes=("_v2", "_full"))
    per_gene = j.groupby("gene").apply(lambda d: pd.Series({
        "r": np.corrcoef(d["mu_log1p_v2"], d["mu_log1p_full"])[0, 1],
        "max_abs": float(np.max(np.abs(d["mu_log1p_v2"] - d["mu_log1p_full"]))),
    }), include_groups=False)
    say(f"  {len(per_gene)} genes; per-gene Pearson r of curves: median {per_gene['r'].median():.4f}, "
        f"5th pct {per_gene['r'].quantile(0.05):.4f}; max |delta log1p|: median {per_gene['max_abs'].median():.4f}, "
        f"95th pct {per_gene['max_abs'].quantile(0.95):.4f}")
    say("\n== Broad/fine modes figure (step 00) vs v1 ==")
    ok = compare_modes()
    (OUT / "comparison_to_v1.txt").write_text("\n".join(lines) + "\n")
    print(f"\nWrote {OUT / 'comparison_to_v1.txt'}")
    if not ok:
        raise SystemExit("modes figure differs from v1")


def compare_modes() -> bool:
    """Step 00 reruns v1's step 01 unchanged on the same r9 tables: its tables must equal v1's exactly."""
    from PIL import Image

    v1m, ok = V1 / "modes", True
    for name in ["dynamic_modes_v3_K562a_fine_cluster_examples_v4.csv", "dynamic_modes_v3_K562a_selected_fine_cluster_examples_v4.csv"]:
        same = pd.read_csv(OUT / name).equals(pd.read_csv(v1m / name))
        ok &= same
        say(f"  {name}: {'identical' if same else 'DIFFERS'}")
    png = "dynamic_modes_v3_K562a_broad_fine_examples_v4.png"
    a, b = (np.asarray(Image.open(p).convert("RGB"), dtype=int) for p in (OUT / png, v1m / png))
    if a.shape == b.shape:
        say(f"  {png}: {a.shape[1]}x{a.shape[0]} px, {int((a != b).any(axis=2).sum())} pixels differ")
    else:
        say(f"  {png}: size differs ({a.shape} vs {b.shape})")
    say("  (v1's check against the original Fig3 figure: input/v1/modes/comparison_to_original.txt)")
    return bool(ok)


if __name__ == "__main__":
    main()
