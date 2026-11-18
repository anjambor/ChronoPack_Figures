"""Axis ranges for the FIG5_UNIVERSE=shared figures (with_prior 5b_v2-5d_v2).

On the 2,077-gene shared universe each TF has fewer reference positives and a higher base rate, so EPR / AUPRCr and
precision values move well outside the full-universe axes (TRRUST AUPRCr up to ~12 vs a 5.5 axis; K562 ChIP values
all near 1 on a 0.9-2.3 axis; PR precision above the 0.35 / 0.5 / 0.7 tops). Only axis limits and ticks change;
values beyond them are still drawn as clipped arrows with their value, as in the full-universe figures.
Called by 02_reproduce_figures.py and 02c_k562_chip_epr_at_500.py when fig5_paths.SHARED_UNIVERSE.
"""
from __future__ import annotations

import numpy as np

# 5b strip-forests: (benchmark, metric) -> (xmin, xmax, ticks)
STRIP_AXES = {
    ("TRRUST", "EPR"): (-0.5, 16.5, [1, 4, 8, 12, 16]),
    ("TRRUST", "AUPRC_ratio"): (0.5, 13.0, [1, 4, 8, 12]),
    ("K562_ChIP", "EPR"): (0.75, 1.7, np.round(np.arange(0.8, 1.6 + 1e-9, 0.2), 1)),
    ("K562_ChIP", "AUPRC_ratio"): (0.85, 1.55, np.round(np.arange(0.9, 1.5 + 1e-9, 0.1), 1)),
}
EPR500_AXIS = (0.85, 1.6, 1.55)  # 02c: xmin, xmax, label threshold
EPR500_TICKS = np.round(np.arange(0.9, 1.5 + 1e-9, 0.1), 1)
# 5c macro-PR: benchmark -> (ylim, yticks); K562 ChIP zoom replaces the full-universe zooms (y 0.2-0.5 / 0.2-0.6)
MACRO_YLIMS = {"TRRUST": ((0.0, 0.4), [0.0, 0.1, 0.2, 0.3, 0.4]), "K562_ChIP": ((0.0, 0.7), [0.0, 0.2, 0.4, 0.6])}
K562_MACRO_ZOOMS = {"_xmax0p5_ymin0p3_ymax0p65": ((0.3, 0.65), [0.3, 0.4, 0.5, 0.6], True)}
# 5d per-TF TRRUST PR curves: y max (the _ymax1p0 variants are unchanged)
PR_YLIMS = {"RELA": 0.8, "NFKB1": 0.8, "JUN": 0.2}


def apply() -> None:
    import plot_publication_v3_6d_trrust_pr_curves as pr6d
    import plot_publication_v3_mean_macro_pr as macro
    import plot_publication_v4_epr_strip_forest as epr_strip

    for (benchmark, metric), (xmin, xmax, ticks) in STRIP_AXES.items():
        metrics = epr_strip.BENCHMARK_CONFIGS[benchmark]["metrics"]
        metrics[:] = [(m, t, xmin, xmax, thr) if m == metric else (m, t, lo, hi, thr) for m, t, lo, hi, thr in metrics]
        epr_strip.X_TICKS[(benchmark, metric)] = ticks
    for benchmark, (ylim, yticks) in MACRO_YLIMS.items():
        macro.YLIMS[benchmark] = ylim
        macro.YTICKS[benchmark] = yticks
    pr6d.YLIMS.update(PR_YLIMS)
