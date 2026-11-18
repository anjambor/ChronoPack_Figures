#!/usr/bin/env python
"""Panel 5d expression curves (RELA, NFKB1, JUN) from the ChronoPack r9 K562 + TNF-alpha
1-lineage tradeSeq fit, as a variant of the capsule's copied (pre-r9, two-curve)
expr-curve_{gene}.pdf panels.

All three genes are refitted together with ChronoPack's quick fit
(cp.tl.fit_gene_dynamics(input_dir=..., genes=[...], final_nknots=8)): the same inputs,
offsets and 8 knots as the r9 final fit, in a temporary folder. RELA is not significant in
r9 (p = 0.65), so it is not in the stored final fit (sce.1trajectory.nknots8.rds); NFKB1 and
JUN are, and their quick-fit curves are checked against it.

Curves: cp.tl.predict_gene_dynamics (fitted mean and 95% band, log1p of the expected count at
the reference library size, ~3,000 UMIs for this run). This is ChronoPack's scale and is lower
than the original panels' (whose normalization is undocumented), so each panel has its own fixed
y-axis (Y_AXES: RELA 0-0.12, ticks every 0.04; NFKB1 and JUN 0-0.3, ticks every 0.1) rather than the
original fixed 0-3 / 0-4. Points: per-time-point mean of normalized counts (cp.tl.observed_gene_dynamics),
drawn as in the original panels. Time tags 1-12 map to 0-110 min (10-min intervals).

Outputs (output/benchmark_outputs/figures/5d/): expr-curve_{gene}_chronopack_1lineage.{pdf,png},
expr-curve_chronopack_1lineage_values.csv, expr-curve_chronopack_1lineage_fit_check.csv.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np
import pandas as pd

import chronopack as cp

from fig5_paths import CHRONOPACK_INPUT_DIR, FIGURE_DIR

FIT_DIR = CHRONOPACK_INPUT_DIR / "K562" / "tradeseq_fit"  # -> r9 sc-K562-1e4_to_tradeSeq.fixed
OUT_DIR = FIGURE_DIR / "5d"
GENES = ["RELA", "NFKB1", "JUN"]
TIME_MAX = 110.0  # K562 TNF-alpha: 12 time tags at 10-min intervals
LABEL = "ChronoPack 1-lineage fit"
COLOR = "#222021"
Y_AXES = {"RELA": (0.12, 0.04), "NFKB1": (0.3, 0.1), "JUN": (0.3, 0.1)}  # gene -> (ymax, tick step); ymin 0
# Measured from the reference PDFs (R/ggplot2, 324 x 324 pt pages): plot area, Helvetica text
# sizes, line widths and mean-point circles, all in points.
REF = dict(panel_w=256.5, panel_h=248.86, tick_pt=20, title_pt=24, tick_pad=9,
           curve_lw=3.63, grid_lw=1.94, mean_edge_lw=1.42, mean_diameter=9.53)


def register_liberation_sans() -> None:
    """Liberation Sans (Helvetica metrics) ships with the env's R package fontLiberation; register it
    so the PNG previews (which cannot use the PDF core Helvetica) render with matching metrics."""
    from matplotlib import font_manager

    font_dir = Path(sys.prefix) / "lib" / "R" / "library" / "fontLiberation" / "fonts" / "liberation-fonts"
    for path in font_dir.glob("LiberationSans-*.ttf"):
        font_manager.fontManager.addfont(str(path))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    register_liberation_sans()
    workers = int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))
    fit = cp.tl.fit_gene_dynamics(input_dir=FIT_DIR, genes=GENES, final_nknots=8, workers=workers, verbose=False)
    curves = cp.tl.predict_gene_dynamics(fit, GENES, n_points=100, time_max=TIME_MAX).assign(lineage=LABEL)
    obs = cp.tl.observed_gene_dynamics(fit, GENES, time_max=TIME_MAX).assign(lineage=LABEL)

    # NFKB1 and JUN are in the stored r9 final fit: the quick fit must reproduce it.
    stored_genes = ["NFKB1", "JUN"]
    stored = cp.tl.predict_gene_dynamics(FIT_DIR, stored_genes, n_points=100, time_max=TIME_MAX,
                                         sce_file="sce.1trajectory.nknots8.rds")
    quick = curves[curves["gene"].isin(stored_genes)].reset_index(drop=True)
    check = pd.DataFrame({
        "gene": stored_genes,
        "max_abs_diff_mu_log1p": [float(np.max(np.abs(quick.loc[quick.gene == g, "mu_log1p"].to_numpy()
                                                      - stored.loc[stored.gene == g, "mu_log1p"].to_numpy())))
                                  for g in stored_genes],
    })
    check.to_csv(OUT_DIR / "expr-curve_chronopack_1lineage_fit_check.csv", index=False)
    print("Quick fit vs stored r9 8-knot fit:\n" + check.to_string(index=False))

    means = (obs.assign(count=np.expm1(obs["expr_log1p"]))
             .groupby(["gene", "t"], as_index=False)["count"].mean()
             .assign(mean_log1p=lambda d: np.log1p(d["count"])).drop(columns="count"))
    values = curves.merge(means, on=["gene", "t"], how="outer").sort_values(["gene", "t"])
    values.to_csv(OUT_DIR / "expr-curve_chronopack_1lineage_values.csv", index=False)

    for gene in GENES:
        # Like the reference PDFs (R pdf device), PDF text uses the core font Helvetica (not embedded;
        # metrics ship with matplotlib). The PNG preview uses Liberation Sans (same metrics).
        with matplotlib.rc_context({"pdf.use14corefonts": True, "font.family": "sans-serif",
                                    "font.sans-serif": ["Helvetica", "Liberation Sans", "DejaVu Sans"]}):
            fig = cp.pl.gene_dynamics_curves(
                curves, [gene], observations=obs, colors={LABEL: COLOR}, ncols=1, panel_size=(6.0, 5.0),
                point_alpha=0.0, y_upper_scale=None, box_aspect=None,
                line_width=REF["curve_lw"], mean_size=REF["mean_diameter"] + REF["mean_edge_lw"],
                mean_edge_width=REF["mean_edge_lw"], tick_label_size=REF["tick_pt"], axis_label_size=REF["title_pt"],
            )
            ax = fig.axes[0]
            for legend in list(fig.legends):
                legend.remove()
            # Same plot-area size (in points) as the reference panels, so text/line/point sizes keep
            # the reference's proportions; the page grows as needed for wider tick labels.
            fig_w, fig_h = fig.get_size_inches()
            ax.set_position([1.6 / fig_w, 0.9 / fig_h, REF["panel_w"] / 72 / fig_w, REF["panel_h"] / 72 / fig_h])
            ax.set_title("", loc="left")  # the original panels have no gene title
            ax.set_xlim(0, TIME_MAX)
            ax.set_xticks([0, 40, 80])
            ymax, step = Y_AXES[gene]
            g = values[values["gene"] == gene]
            if float(np.nanmax(g[["hi_log1p", "mean_log1p"]].to_numpy())) > ymax:
                raise ValueError(f"{gene}: curve band or means exceed ymax {ymax}")
            ax.set_ylim(0, ymax)
            ax.set_yticks(np.round(np.arange(0, ymax + step / 2, step), 10))
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:g}"))
            ax.grid(True, which="major", color="#e0e0e0", lw=REF["grid_lw"], zorder=0)
            ax.tick_params(labelsize=REF["tick_pt"], pad=REF["tick_pad"])
            ax.set_xlabel("Time (min)", fontsize=REF["title_pt"], labelpad=4)
            ax.set_ylabel("Expression", fontsize=REF["title_pt"], labelpad=8)
            stem = OUT_DIR / f"expr-curve_{gene}_chronopack_1lineage"
            fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.1)
            fig.savefig(stem.with_suffix(".png"), bbox_inches="tight", pad_inches=0.1, dpi=300)
            plt.close(fig)
        print(f"Wrote {stem}.pdf/.png")

if __name__ == "__main__":
    main()
