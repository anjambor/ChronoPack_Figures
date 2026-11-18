#!/usr/bin/env python
"""Figures for the P < 0.05 BEELINE round (from step 10's output/rounds/per_tf_metrics.csv).

1. p05_strip_forest_<benchmark>_full_TFs<N>.pdf: Fig5 panel 5b (per-TF EPR and AUPRCr strip-forest, full Fig5
   universe) with the K562 BEELINE rows replaced by the r9 P < 0.05 runs and THP1 BEELINE rows added (THP1 PPCOR is
   not applicable). ChronoPack, Replogle and Yao rows are Fig5's.
2. p05_strip_forest_<benchmark>_matched_TFs<N>.pdf: the same for 04e's matched-gene universe (each context's BEELINE
   genes also ranked by that context's ChronoPack), ChronoPack vs BEELINE only.
3. rounds_median_<metric>_TFs<N>.pdf: median over the 5b TFs per round (10% filter, no filter, P < 0.05; the metacell
   round of Fig5 v1 was dropped 2026-11-16),
   full and matched universes, both contexts and benchmarks: how gene selection moves the BEELINE scores.

Drawing uses Fig5's 5b code from Fig5_GRN (v2's, copied unchanged through v3; plot_strip_forest_v3: draw_metric,
spread_dots, place_labels), imported read-only; only the method rows, labels and axis ranges are set here. Writes output/rounds/figures/.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

import os  # noqa: E402
os.environ["FIG5_METHOD_SET"] = "expression"
V2 = c.ROOT.parent / "scripts"  # Fig5_GRN/scripts: 5b plotting code (v2's, copied unchanged through v3 into v4)
for p in (V2, V2 / "capsule_lib/01_generate_predictions_and_metrics", V2 / "capsule_lib/02_visualize_publication_figures"):
    sys.path.insert(0, str(p))
OUT = c.OUTPUT_DIR / "rounds" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".mplconfig"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import benchmark_grn_systems as b  # noqa: E402
import plot_publication_v3_trrust_strip_forest as strip  # noqa: E402
import plot_publication_v4_epr_strip_forest as v2  # noqa: E402
import plot_strip_forest_v3 as v3  # noqa: E402

BEELINE_K562 = ["GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"]
BEELINE_THP1 = ["GENIE3_THP1", "SINCERITIES_THP1", "GRNBOOST2_THP1"]  # PPCOR not applicable for THP1
for m in BEELINE_THP1:
    base = m.removesuffix("_THP1")
    b.METHOD_COLORS[m] = b.METHOD_COLORS[base]
    strip.VIOLIN_METHOD_LABELS[m] = f"{base}\nTHP1"
v2.THP1_METHODS = set(v2.THP1_METHODS) | set(BEELINE_THP1)
strip.ordered_tfs_for_benchmark = v2.ordered_tfs_for_benchmark
FULL_ROWS = ["ChronoPack_K562", *BEELINE_K562, "Replogle2022_K562", "ChronoPack_THP1", *BEELINE_THP1, "Yao2023_THP1"]
MATCHED_ROWS = ["ChronoPack_K562", *BEELINE_K562, "ChronoPack_THP1", *BEELINE_THP1]
# x ranges (xmin, xmax) per (figure, benchmark, metric); full = Fig5 5b's ranges (v2.BENCHMARK_CONFIGS)
RANGES = {
    ("matched", "TRRUST", "EPR"): (-1.0, 25.5), ("matched", "TRRUST", "AUPRC_ratio"): (0.0, 18.0),
    ("matched", "K562_ChIP", "EPR"): (0.6, 1.6), ("matched", "K562_ChIP", "AUPRC_ratio"): (0.6, 1.6),
}
TITLE = {"TRRUST": "TRRUST", "K562_ChIP": "K562 ChIP"}


def strip_forest(data: pd.DataFrame, rows: list[str], benchmark: str, kind: str, n: int) -> Path:
    b.METHOD_ORDER = rows
    strip.METHOD_ORDER = rows
    data = v2.filter_rows(data[data["benchmark"].eq(benchmark) & data["method"].isin(rows)].copy())
    gap = 0.5  # extra space between the K562 and THP1 blocks
    first_thp1 = rows.index("ChronoPack_THP1")
    y = {m: 0.86 * i + (gap if i >= first_thp1 else 0) for i, m in enumerate(rows)}
    height = {"TRRUST": 6.8, "K562_ChIP": 4.65}[benchmark] * (len(rows) / 8) ** 0.85
    fig, axes = plt.subplots(1, 2, figsize=(10.6, height), sharey=True)
    clipped, summary, drawn = [], [], []
    for ax, (metric, title, xmin, xmax, _t) in zip(axes, v2.BENCHMARK_CONFIGS[benchmark]["metrics"]):
        xmin, xmax = RANGES.get((kind, benchmark, metric), (xmin, xmax))
        drawn.append(v3.draw_metric(ax, data, benchmark, metric, title, xmin, xmax, y, clipped, summary))
        if kind == "matched":
            ax.set_xticks(np.arange(0, xmax + 1e-9, 5 if benchmark == "TRRUST" else 0.2))
            if benchmark == "K562_ChIP":
                ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))
    for ax in axes:
        ax.set_yticks([y[m] for m in rows])
        ax.set_yticklabels([strip.method_label(m) for m in rows], fontsize=11.0)
        ax.set_ylim(max(y.values()) + 0.55, -0.62)
    ymax = max(y.values()) + 0.55 + 0.62
    for label, members in (("K562", rows[:first_thp1]), ("THP1", rows[first_thp1:])):
        mid = (y[members[0]] + y[members[-1]]) / 2
        axes[0].text(-0.30, 1 - (mid + 0.62) / ymax, label, transform=axes[0].transAxes, ha="center", va="center",
                     rotation=90, fontsize=9.8, color="#555555")
    legend = axes[1].legend(handles=[Line2D([0], [0], marker="D", linestyle="", markerfacecolor="white",
                                            markeredgecolor="#111111", markersize=7.8, label="Median")],
                            frameon=False, loc="lower right", fontsize=10.0, handletextpad=0.4)
    universe = "Fig. 5 universe" if kind == "full" else "matched genes (BEELINE genes ranked by ChronoPack)"
    fig.suptitle(f"{TITLE[benchmark]}: BEELINE P < 0.05, TFs + {n}; {universe}", fontsize=14.0, y=0.985)
    fig.tight_layout(rect=[0.055, 0.04, 1, 0.945], w_pad=2.0)
    fig.text(0.01, 0.005, "THP1 PPCOR: not applicable (787 cells < genes).", fontsize=8.0, color="#555555")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    warnings = []
    for ax, (points, medians) in zip(axes, drawn):
        v3.spread_dots(fig, ax, points, medians, *v3.DOT_SPREAD[benchmark])
    for i, (ax, (points, medians)) in enumerate(zip(axes, drawn)):
        extra = [tuple(p["annotation"].get_window_extent(renderer).extents) for p in points if "annotation" in p]
        if i == 1:
            extra.append(tuple(legend.get_window_extent(renderer).extents))
        warnings += v3.place_labels(fig, ax, points, medians, extra, refine=v3.REFINE_LABELS[benchmark],
                                    row_band=v3.ROW_BAND[benchmark])
    out = OUT / f"p05_strip_forest_{benchmark}_{kind}_TFs{n}.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(summary).to_csv(out.with_suffix(".medians.csv"), index=False)
    pd.DataFrame(clipped).to_csv(out.with_suffix(".clipped.csv"), index=False)
    print(out, *(f"\n  unresolved: {w}" for w in warnings))
    return out


ROUND_ORDER = ["r9_10pct", "r9_nofilter", "r9_p05"]  # metacell_Fig5 dropped 2026-11-16
ROUND_LABEL = {"metacell_Fig5": "metacell\n(Fig5 v1)", "r9_10pct": "10%\nfilter", "r9_nofilter": "no\nfilter",
               "r9_p05": "P < 0.05"}


def rounds_figure(med: pd.DataFrame, metric: str, n: int) -> Path:
    col = {"EPR": "median_EPR", "AUPRCr": "median_AUPRCr"}[metric]
    panels = [(ctx, bench) for ctx in ("K562", "THP1") for bench in ("TRRUST", "K562_ChIP")]
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.2), sharex=False)
    for j, (ctx, bench) in enumerate(panels):
        rounds = [r for r in ROUND_ORDER if ctx == "K562" or r != "metacell_Fig5"]
        xs = np.arange(len(rounds))
        methods = [f"ChronoPack_{ctx}"] + (BEELINE_K562 if ctx == "K562" else BEELINE_THP1)
        for i, uni in enumerate(["full", "shared"]):
            ax = axes[i, j]
            sub = med[(med.context == ctx) & (med.benchmark == bench) & (med.universe == uni)
                      & (med.TFsN.isna() | med.TFsN.eq(n))]
            for m in methods:
                v = [sub.loc[(sub["round"] == r) & (sub.method == m), col] for r in rounds]
                v = [float(s.iloc[0]) if len(s) else np.nan for s in v]
                cp = m.startswith("ChronoPack")
                ax.plot(xs, v, marker="D" if cp else "o", ms=7 if cp else 5.5, lw=2.2 if cp else 1.4,
                        color=b.METHOD_COLORS[m], label=strip.method_label(m).split(" |")[0], zorder=4 if cp else 3)
            ax.axhline(1, color="#9a9a9a", lw=0.9, ls="--", zorder=0)
            ax.set_xticks(xs)
            ax.set_xticklabels([ROUND_LABEL[r] for r in rounds], fontsize=8.5)
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="y", color="#e7e7e7", lw=0.6)
            if i == 0:
                ax.set_title(f"{ctx} · {TITLE[bench]}", fontsize=11)
            if j == 0:
                ax.set_ylabel(f"median {metric}\n{'Fig. 5 universe' if uni == 'full' else 'matched genes'}", fontsize=10)
    for j, x in ((0, 0.26), (2, 0.76)):  # one legend per context, below its two columns
        h, l = axes[0, j].get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", bbox_to_anchor=(x, 0.035), ncol=len(h), fontsize=9, frameon=False)
    fig.suptitle(f"BEELINE gene selection rounds, TFs + {n}: median {metric} over the 5b TFs "
                 f"(ChronoPack is the same predictions in every round)", fontsize=12.5)
    have = med[med.universe.eq("full") & (med.TFsN.isna() | med.TFsN.eq(n))]
    missing = [f"{ctx} {m.removesuffix('_THP1')} ({ROUND_LABEL[r].replace(chr(10), ' ')})"
               for ctx, bl in (("K562", BEELINE_K562), ("THP1", BEELINE_THP1)) for r in ROUND_ORDER
               if not (ctx == "THP1" and r == "metacell_Fig5") for m in bl
               if have[(have.context == ctx) & (have["round"] == r) & (have.method == m)].empty]
    fig.text(0.01, 0.005, "THP1 PPCOR: not applicable (787 cells < genes). Matched genes: each round's BEELINE genes "
             "also ranked by ChronoPack of that context."
             + (f" Not finished (gap): {', '.join(missing)}." if missing else ""), fontsize=8.0, color="#555555")
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    out = OUT / f"rounds_median_{metric}_TFs{n}.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(out)
    return out


def main() -> None:
    m = pd.read_csv(c.OUTPUT_DIR / "rounds" / "per_tf_metrics.csv")
    m["EPR"] = m["EPR_capped"]  # 5b plots the capped EPR
    med = pd.read_csv(c.OUTPUT_DIR / "rounds" / "median_by_universe.csv")
    for n in c.N_GENES:
        sel = m[m["round"].eq("r9_p05") & m["TFsN"].eq(n)]
        others = sel.method.isin(["ChronoPack_K562", "Replogle2022_K562", "ChronoPack_THP1", "Yao2023_THP1"])
        # ChronoPack / Replogle / Yao: one copy (identical in every run of a context in the full universe)
        full = pd.concat([sel[sel.universe.eq("full") & ~others],
                          sel[sel.universe.eq("full") & others].drop_duplicates(["benchmark", "method", "TF"])])
        matched = sel[sel.universe.eq("shared")]
        for bench in ("TRRUST", "K562_ChIP"):
            strip_forest(full, FULL_ROWS, bench, "full", n)
            strip_forest(matched, MATCHED_ROWS, bench, "matched", n)
        for metric in ("AUPRCr", "EPR"):
            rounds_figure(med, metric, n)
    print("P05_FIGURES_OK")


if __name__ == "__main__":
    main()
