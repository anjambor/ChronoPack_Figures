#!/usr/bin/env python
"""Fig5 composite (panels a-e) laid out as ref/FigureComposite.v6.png, built from Fig5_GRN_v4's expression method set.

ref/FigureComposite.v6.png is the user's composite of the metacell-era Fig5 (5,000 x 4,111 px, 300 dpi). This step
redraws its subfigures from v3's tables and places each one in the bounding box the reference panel occupies
(REF_BOXES, measured from the reference's non-white pixels; panel (a), the method schematic, is left blank):

  (b) 5b strip-forests, TRRUST and K562 ChIP, with the reference's styling: y-axis lines, bracketed K562 / THP-1
      groups, bold title and axis labels, Median legend in ChIP only, and sparse TF labels beside their dots
      (strip_label_layout; COMPOSITE_STRIP_LABELS=all for v3's every-TF placement).
  (c) 5c mean macro-PR curves: TRRUST (x 0-0.25) with the legend, and the K562 ChIP zoom (x 0-0.5, y 0.2-0.6)
      with its title inside the axes.
  (d) 5e rank-enrichment heatmap, ChronoPack and the BEELINE methods (no fold column), colour bar beside it.
  (e) 5d TRRUST heatmap (RELA / NFKB1 / JUN), the per-TF TRRUST PR curves (y up to 1.0 for RELA and NFKB1, 0.1 for
      JUN; no legend) and the ChronoPack expression curves (02b's PNG / PDF, unchanged).

All text is Liberation Sans (Arial's metrics), as in the reference. Each panel is rendered at 600 dpi, cropped to its
non-white pixels (as REF_BOXES were measured), its figure size adjusted until that area has the reference panel's
aspect ratio, then scaled into the reference box. Data are the expression method set's tables (FIG5_METHOD_SET=
expression; run 02 and 02b first). The PNG is the raster canvas (5,000 x 4,111 px, 300 dpi); the PDF is vector: each
panel's PDF page, cropped the same way, placed in the same box (vector_pdf; needs pypdf in Fig5_GRN/.pylib:
`python -m pip install --no-deps --target .pylib pypdf`). Text in the PDF is embedded TrueType.
Writes output/composite/panels/*.{png,pdf} and output/composite/Fig5_composite_v4.{png,pdf}.
"""
from __future__ import annotations

import io
import os
import sys
from pathlib import Path

os.environ["FIG5_METHOD_SET"] = "expression"
SCRIPTS = Path(__file__).resolve().parent
for p in (SCRIPTS, SCRIPTS / "capsule_lib" / "01_generate_predictions_and_metrics", SCRIPTS / "capsule_lib" / "02_visualize_publication_figures"):
    sys.path.insert(0, str(p))

import fig5_paths  # noqa: E402

os.environ.setdefault("MPLCONFIGDIR", str(fig5_paths.OUTPUT_DIR / ".mplconfig"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
from matplotlib import font_manager  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.transforms as mtransforms  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.ticker import FormatStrFormatter  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402

FONT = "Liberation Sans"
# Liberation Sans (Arial's metrics) ships with the env's R package fontLiberation; register it as 02b does (compute
# nodes have no system copy)
for _f in (Path(sys.prefix) / "lib" / "R" / "library" / "fontLiberation" / "fonts" / "liberation-fonts").glob("LiberationSans-*.ttf"):
    font_manager.fontManager.addfont(str(_f))
_set_theme = sns.set_theme
sns.set_theme = lambda *a, **k: _set_theme(*a, **{**k, "font": FONT})  # every module's theme call keeps the font

import benchmark_grn_systems as b  # noqa: E402
import generate_system_broad_filtered_publication_figures as pubfig  # noqa: E402
import plot_publication_v3_mean_macro_pr as macro  # noqa: E402
import plot_publication_v3_trrust_core_heatmap_reference_palette as heat  # noqa: E402
import plot_publication_v3_trrust_strip_forest as strip  # noqa: E402
import plot_publication_v4_epr_strip_forest as v2  # noqa: E402
import plot_regulator_prioritization_v2 as reg  # noqa: E402
import plot_regulator_prioritization_v2_composite_stacked as comp  # noqa: E402
import plot_regulator_prioritization_v2_nes_dot_comparison as nes  # noqa: E402
import plot_strip_forest_v3 as v3  # noqa: E402
import strip_label_layout as sparse  # noqa: E402
from generate_publication_figures_v4_filtered_tfs import K562_DROP_TFS, filter_curve_rows  # noqa: E402

OUT = fig5_paths.OUTPUT_DIR / "composite"
PANELS = OUT / "panels"
REF = fig5_paths.FIG_DIR / "ref" / "FigureComposite.v6.png"
FIG5D = fig5_paths.FIGURE_DIR / "5d"
CANVAS = (5000, 4111)
DPI_OUT = 300
RENDER_DPI = 600
# Content bounding boxes (x0, y0, x1, y1, px) of the reference panels, from its non-white pixels
REF_BOXES = {
    "b_trrust": (2259, 56, 4818, 1150), "b_chip": (2241, 1212, 4846, 2305),
    "c_trrust": (2380, 2414, 3194, 3224), "c_chip": (2376, 3328, 3187, 4073),
    "d": (3440, 2397, 4812, 3995),
    "e_heat": (266, 2113, 2076, 2807),
    "e_pr_RELA": (272, 2890, 790, 3461), "e_pr_NFKB1": (929, 2890, 1447, 3461), "e_pr_JUN": (1577, 2885, 2114, 3461),
    "e_expr_RELA": (273, 3522, 803, 4021), "e_expr_NFKB1": (938, 3522, 1450, 4021), "e_expr_JUN": (1594, 3522, 2106, 4021),
}
LETTERS = {"(a)": (68, 50), "(b)": (2194, 50), "(c)": (2194, 2361), "(d)": (3344, 2361), "(e)": (68, 2046)}  # top-left
LETTER_HEIGHT = 95  # px, "(" top to ")" bottom in the reference
STRIP_TOP = {"TRRUST": 0.978, "K562_ChIP": 0.993}  # strip-forest axes top (figure fraction): reference title gaps
STRIP_XMIN = {("TRRUST", "EPR"): -0.75}  # reference: EPR = 0 dots farther from the y axis than 5b's -0.5
GROUP_COLOR = "#2d4056"
THP1_LABEL = "THP-1"
# (b) TF labels: "sparse" = the reference's style (labels beside their dots, dots nudged within rows, labels left out
# where none fits; strip_label_layout); "all" = every TF labelled with v3's placement (5b's standalone figures)
STRIP_LABELS = os.environ.get("COMPOSITE_STRIP_LABELS", "sparse")
_SIZING = False  # True while render_to_aspect searches the figure size (labels skipped: they stay inside the axes)


def style() -> None:
    """Reference look: seaborn white theme, Liberation Sans, black text and ticks."""
    sns.set_theme(style="white", font_scale=0.9)
    plt.rcParams.update({"pdf.fonttype": 42, "font.family": FONT, "text.color": "black", "axes.labelcolor": "black",
                         "xtick.color": "black", "ytick.color": "black", "axes.edgecolor": "black"})


# ---------------------------------------------------------------- rendering helpers
def tight_png(fig, dpi: int = RENDER_DPI) -> Image.Image:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", pad_inches=0.02, facecolor="white")
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def ink_box(im: Image.Image) -> tuple[int, int, int, int]:
    """Bounding box of the non-white pixels (REF_BOXES are measured the same way)."""
    a = np.asarray(im.convert("L")) < 250
    ys, xs = np.where(a.any(1))[0], np.where(a.any(0))[0]
    return int(xs[0]), int(ys[0]), int(xs[-1]) + 1, int(ys[-1]) + 1


def as_panel(full: Image.Image, pdf: Path) -> dict:
    """Panel cropped to its ink; frac = the crop as fractions of the full image (= of the PDF page, same bbox)."""
    box = ink_box(full)
    frac = (box[0] / full.width, box[1] / full.height, box[2] / full.width, box[3] / full.height)
    return {"im": full.crop(box), "pdf": pdf, "frac": frac}


def render_to_aspect(draw, figsize: tuple[float, float], aspect: float, stem: str, iters: int = 6) -> dict:
    """draw(figsize) -> Figure. Adjusts the figure width until the inked area has the target aspect."""
    global _SIZING
    w, h = figsize
    _SIZING = True
    for _ in range(iters):
        fig = draw((w, h))
        im = tight_png(fig)  # at the final resolution: the ink edge (anti-aliasing) depends on it
        plt.close(fig)
        im = im.crop(ink_box(im))
        got = im.width / im.height
        if abs(got / aspect - 1) < 0.002:
            break
        w += (aspect * im.height - im.width) / RENDER_DPI
    _SIZING = False
    fig = draw((w, h))
    # transparent background: the PDF's padding must not paint over neighbouring panels in the vector composite
    fig.savefig(PANELS / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.02, transparent=True)
    panel = as_panel(tight_png(fig), PANELS / f"{stem}.pdf")
    plt.close(fig)
    panel["im"].save(PANELS / f"{stem}.png")
    print(f"  {stem}: figsize {w:.2f} x {h:.2f} in, ink aspect {panel['im'].width / panel['im'].height:.3f} "
          f"(target {aspect:.3f})", flush=True)
    return panel


def box_aspect(key: str) -> float:
    x0, y0, x1, y1 = REF_BOXES[key]
    return (x1 - x0) / (y1 - y0)


# ---------------------------------------------------------------- (b) strip-forests
def strip_panel(data: pd.DataFrame, benchmark: str, figsize) -> plt.Figure:
    style()
    v3.LABEL_FONTSIZE = 8.6  # reference: TF labels ~13% larger than v3's 7.6
    rows = b.METHOD_ORDER
    sub = data[data["benchmark"].eq(benchmark) & data["method"].isin(rows)].copy()
    y_values, _h, _labels = b.strip_rows(rows, figsize[1])
    y = dict(zip(rows, y_values))
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    clipped, summary, drawn = [], [], []
    for ax, (metric, title, xmin, xmax, _t) in zip(axes, v2.BENCHMARK_CONFIGS[benchmark]["metrics"]):
        xmin = STRIP_XMIN.get((benchmark, metric), xmin)
        drawn.append(v3.draw_metric(ax, sub, benchmark, metric, title, xmin, xmax, y, clipped, summary))
        ax.set_xlabel(title, fontsize=13.5, fontweight="bold")
        ax.tick_params(axis="x", labelsize=12, length=4.0, pad=2.2)  # reference tick length and label gap
        ax.spines["left"].set_visible(True)  # reference: y-axis line on both panels
        ax.spines["left"].set_color(ax.spines["bottom"].get_edgecolor())
        ax.spines["left"].set_linewidth(ax.spines["bottom"].get_linewidth())
    for ax in axes:
        ax.set_yticks([y[m] for m in rows])
        ax.set_yticklabels([strip.method_label(m) for m in rows], fontsize=12.5)
        ax.set_ylim(max(y_values) + 0.55, -0.62)
    legend = None
    if benchmark == "K562_ChIP":  # reference: Median legend in the ChIP panel only
        legend = axes[1].legend(handles=[Line2D([0], [0], marker="D", linestyle="", markerfacecolor="white",
                                                markeredgecolor="#111111", markersize=7.8, label="Median")],
                                frameon=False, loc="lower right", fontsize=10.0, handletextpad=0.4)
    fig.suptitle(v2.BENCHMARK_CONFIGS[benchmark]["label"], fontsize=24, fontweight="bold", y=1.0, va="top")
    fig.tight_layout(rect=[0.075, 0.02, 1, STRIP_TOP[benchmark]], w_pad=2.0)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    if STRIP_LABELS == "all":
        for ax, (points, medians) in zip(axes, drawn):
            v3.spread_dots(fig, ax, points, medians, *v3.DOT_SPREAD[benchmark])
    records = []
    for i, (ax, (points, medians)) in enumerate(zip(axes, drawn)):
        if _SIZING:
            break
        extra = [tuple(p["annotation"].get_window_extent(renderer).extents) for p in points if "annotation" in p]
        if i == 1 and legend is not None:
            extra.append(tuple(legend.get_window_extent(renderer).extents))
        if STRIP_LABELS == "all":
            v3.place_labels(fig, ax, points, medians, extra, refine=v3.REFINE_LABELS[benchmark], row_band=v3.ROW_BAND[benchmark])
        else:
            metric = v2.BENCHMARK_CONFIGS[benchmark]["metrics"][i][1]
            records += [{"benchmark": benchmark, "metric": metric, **r} for r in
                        sparse.layout(fig, ax, points, medians, extra, fontsize=v3.LABEL_FONTSIZE, family=FONT)]
    if records:
        rec = pd.DataFrame(records)
        rec.to_csv(PANELS / f"b_{'trrust' if benchmark == 'TRRUST' else 'chip'}_labels.csv", index=False)
        print(f"  {benchmark}: {int(rec['labelled'].astype(bool).sum())} of {len(rec)} dots labelled; "
              f"max residual cost {rec['residual_cost'].max():.1f}", flush=True)
    add_strip_groups(fig, axes[0], rows, y)
    return fig


def add_strip_groups(fig, ax, rows, y) -> None:
    """Reference style: a vertical line left of each block's method labels, the group name rotated beside it."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = ax.transAxes.inverted()
    ticks = {lab.get_text(): lab for lab in ax.get_yticklabels()}
    labels = ax.get_yticklabels()
    blocks = [("K562", [m for m in rows if m not in b.THP1_METHODS]), (THP1_LABEL, [m for m in rows if m in b.THP1_METHODS])]
    for name, members in blocks:
        idx = [rows.index(m) for m in members]
        left = min(labels[i].get_window_extent(renderer).x0 for i in idx)
        x_line = inv.transform((left - 9, 0))[0]
        trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
        y0, y1 = y[members[0]] - 0.38, y[members[-1]] + 0.38
        ax.plot([x_line, x_line], [y0, y1], transform=trans, color=GROUP_COLOR, lw=1.6, clip_on=False, solid_capstyle="butt")
        x_text = inv.transform((left - 22, 0))[0]
        ax.text(x_text, (y0 + y1) / 2, name, transform=trans, rotation=90, ha="right", va="center", fontsize=13,
                fontweight="bold", color=GROUP_COLOR, clip_on=False)


# ---------------------------------------------------------------- (c) macro-PR
MACRO_LEGEND = ["ChronoPack_THP1", "ChronoPack_K562", "GENIE3", "Yao2023_THP1", "Replogle2022_K562"]
MACRO_LABELS = {"ChronoPack_THP1": "Chrono | THP1", "ChronoPack_K562": "Chrono | K562", "GENIE3": "GENIE3 | K562",
                "Yao2023_THP1": "Yao | THP1", "Replogle2022_K562": "Replogle | K562"}


def mean_curve(curves: pd.DataFrame, benchmark: str, method: str):
    """macro.plot_macro_pr_for_benchmark's mean precision over TFs on a 250-point recall grid (unchanged logic)."""
    grid = np.linspace(0, 1, 250)
    panel = curves[(curves["benchmark"].eq(benchmark)) & (curves["method"].eq(method))]
    ys, masks = [], []
    for tf in sorted(panel["TF"].dropna().unique()):
        d = panel[panel["TF"].eq(tf)].sort_values(["recall"]).dropna(subset=["recall", "precision"]).drop_duplicates("recall", keep="last")
        x, yv = d["recall"].to_numpy(float), d["precision"].to_numpy(float)
        if len(x) == 0:
            continue
        o = np.argsort(x)
        ys.append(np.interp(grid, x[o], yv[o], left=yv[o][0], right=yv[o][-1]))
        obs = d[~d["is_expectation"]]
        masks.append(grid > (obs["recall"].max() if len(obs) else -np.inf))
    if not ys:
        return None
    return grid, np.nanmean(np.vstack(ys), axis=0), np.any(np.vstack(masks), axis=0)


def segmented(ax, x, y, dashed, color, label) -> None:
    """macro.plot_segmented_line with the reference's thicker lines (solid 3.4, dashed 2.8)."""
    used = False
    for target, ls, lw, alpha in [(False, "-", 3.4, 1.0), (True, "--", 2.8, 0.72)]:
        idx = np.flatnonzero(dashed == target)
        if not len(idx):
            continue
        for g in np.split(idx, np.where(np.diff(idx) > 1)[0] + 1):
            if len(g) == 1:
                g = np.array([max(0, g[0] - 1), g[0]])
            ax.plot(x[g], y[g], lw=lw, ls=ls, color=color, alpha=alpha, label=None if used else label)
            used = True


def macro_panel(curves: pd.DataFrame, benchmark: str, figsize) -> plt.Figure:
    style()
    fig, ax = plt.subplots(figsize=figsize)
    for method in pubfig.PR_METHOD_ORDER:
        r = mean_curve(curves, benchmark, method)
        if r is not None:
            segmented(ax, *r, color=b.METHOD_COLORS[method], label=MACRO_LABELS[method])
    if benchmark == "TRRUST":
        ax.set_xlim(0, 0.25)
        ax.set_ylim(0, 0.35)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3])
        ax.set_title("TRRUST", fontsize=24, fontweight="bold", pad=10)
        h, l = ax.get_legend_handles_labels()
        order = [l.index(MACRO_LABELS[m]) for m in MACRO_LEGEND if MACRO_LABELS[m] in l]
        ax.legend([h[i] for i in order], [l[i] for i in order], frameon=False, fontsize=11.5, loc="upper right",
                  handlelength=1.6, borderaxespad=0.6, labelspacing=0.35)
    else:
        ax.set_xlim(0, 0.5)
        ax.set_ylim(0.2, 0.6)
        ax.set_xticks([0.0, 0.5])
        ax.set_yticks([0.2, 0.4, 0.6])
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.text(0.5, 0.98, "K562 ChIP", transform=ax.transAxes, ha="center", va="top", fontsize=24, fontweight="bold")
    ax.set_xlabel("Recall", fontsize=18)
    ax.set_ylabel("Precision", fontsize=18)
    ax.tick_params(axis="both", labelsize=18)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- (d) rank-enrichment heatmap
def nes_panel(figsize) -> plt.Figure:
    style()
    nes.K562_TFS = [tf for tf in reg.K562_EVAL_TFS if tf not in K562_DROP_TFS]  # as v4.generate_6e (5e)
    matrix, summary = comp.load_k562_matrix()
    log_values = np.log10(matrix.to_numpy(float))
    finite = log_values[np.isfinite(log_values)]
    norm = Normalize(vmin=float(finite.min()), vmax=float(finite.max()))
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.15, 0.055], wspace=0.0, left=0.1, right=0.86, top=0.86, bottom=0.03)
    ax, cax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 2])
    n_rows, n_cols = matrix.shape
    # pcolormesh, not imshow: cells stay vector shapes in the PDF
    im = ax.pcolormesh(np.arange(n_cols + 1) - 0.5, np.arange(n_rows + 1) - 0.5, np.ma.masked_invalid(log_values),
                       cmap=comp.MILD_VIRIDIS, norm=norm, shading="flat", linewidth=0)
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([comp.NES_METHOD_LABELS[m] for m in matrix.columns], rotation=40, ha="left", va="bottom",
                       rotation_mode="anchor", fontsize=12.5)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(matrix.index.tolist(), fontsize=13.5)
    ax.tick_params(axis="both", length=0, pad=4)
    ax.xaxis.tick_top()
    for i, tf in enumerate(matrix.index):
        for j, method in enumerate(matrix.columns):
            v = matrix.loc[tf, method]
            ax.text(j, i, comp.format_nes(v), ha="center", va="center", fontsize=13,
                    color=comp.text_color(np.log10(v), norm), fontweight="bold" if method == "ChronoPack_K562" else "normal")
    for x in np.arange(-0.5, n_cols, 1):
        ax.axvline(x, color="white", lw=1.2)
    for yy in np.arange(-0.5, n_rows, 1):
        ax.axhline(yy, color="white", lw=1.2)
    ax.add_patch(Rectangle((-0.5, -0.5), 1, n_rows, fill=False, edgecolor="#2b2b2b", linewidth=1.6, zorder=5))
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.spines[:].set_visible(False)
    # colour bar: beside the heatmap, spanning the middle ~70% of it, ticks and label on the right (reference)
    pos = ax.get_position()
    cpos = cax.get_position()
    cax.set_position([cpos.x0, pos.y0 + pos.height * 0.07, cpos.width, pos.height * 0.86])
    cbar = plt.colorbar(im, cax=cax, orientation="vertical")
    cbar.solids.set_rasterized(False)  # matplotlib rasterizes colour bars with >= 50 steps; keep it vector
    ticks = np.array([0.01, 0.1, 1, 10, 100])
    ticks = ticks[(np.log10(ticks) >= norm.vmin) & (np.log10(ticks) <= norm.vmax)]
    cbar.set_ticks(np.log10(ticks))
    cbar.set_ticklabels([comp.log_tick_formatter(t, 0) for t in ticks])
    cbar.outline.set_linewidth(1.6)
    cbar.outline.set_edgecolor("#2b2b2b")
    cbar.ax.tick_params(labelsize=17, length=4, pad=3)
    cbar.set_label("Rank enrichment score (log scale)", fontsize=17, rotation=270, labelpad=24)
    return fig


# ---------------------------------------------------------------- (e) heatmap and per-TF PR
def heat_panel(figsize) -> plt.Figure:
    """heat.main's figure (5d TRRUST RELA / NFKB1 / JUN heatmap) with the reference's THP-1 label, returned unsaved."""
    captured = {}

    def keep(fig):
        captured["fig"] = fig

    def groups(ax) -> None:  # pubfig.add_trrust_tf_cell_type_row_groups with "THP-1"
        trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
        n_k562 = len([m for m in pubfig.PR_METHOD_ORDER if m not in b.THP1_METHODS])
        for label, y0, y1 in [("K562", 0, n_k562), (THP1_LABEL, n_k562, len(pubfig.PR_METHOD_ORDER))]:
            ax.text(-0.27, (y0 + y1) / 2, label, transform=trans, ha="center", va="center", rotation=90, fontsize=11.5,
                    fontweight="bold", color=GROUP_COLOR, clip_on=False)
            ax.plot([-0.235, -0.215, -0.215, -0.235], [y0 + 0.08, y0 + 0.08, y1 - 0.08, y1 - 0.08], transform=trans,
                    color=GROUP_COLOR, lw=1.4, clip_on=False)

    heat.FIGSIZE = figsize
    heat.SOURCE_OUT = PANELS / "e_heat_source.csv"
    heat.add_trrust_tf_cell_type_row_groups = groups
    saved_save = heat.save
    heat.save = keep
    try:
        heat.main()
    finally:
        heat.save = saved_save
    return captured["fig"]


def pr_panel(curves: pd.DataFrame, tf: str, figsize) -> plt.Figure:
    import plot_publication_v3_6d_trrust_pr_curves as pr6d
    style()
    panel = curves[curves["scenario"].eq(pr6d.SCENARIO) & curves["benchmark"].eq("TRRUST") & curves["scope"].eq("per_tf")
                   & curves["TF"].eq(tf) & curves["method"].isin(pr6d.METHOD_ORDER)]
    ymax = pr6d.VARIANT_YLIM if tf in ("RELA", "NFKB1") else pr6d.YLIMS[tf]
    fig, ax = plt.subplots(figsize=figsize)
    for method in pr6d.METHOD_ORDER:
        d = panel[panel["method"].eq(method)].sort_values("k")
        obs = d[~d["is_expectation"]].dropna(subset=["recall", "precision"])
        exp = d[d["is_expectation"]].dropna(subset=["recall", "precision"])
        color = pr6d.METHOD_COLORS[method]
        if len(obs):
            ax.plot(obs["recall"], obs["precision"], drawstyle="steps-post", lw=2.0, color=color)
        if len(exp):
            start = obs.tail(1) if len(obs) else exp.head(1)
            tail = pd.concat([start, exp]).drop_duplicates(["recall", "precision"])
            ax.plot(tail["recall"], tail["precision"], lw=1.7, ls="--", color=color, alpha=0.75)
    ax.set_xlim(0, 0.25)
    ax.set_ylim(0, ymax)
    if tf == "JUN":
        ax.set_yticks(np.arange(0, ymax + 0.001, 0.04))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    else:
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.set_xticks([0.0, 0.1, 0.2])
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.set_box_aspect(1)
    ax.set_xlabel("Recall", fontsize=18, labelpad=2)
    ax.set_ylabel("Precision", fontsize=18, labelpad=3)
    ax.set_title(tf, fontsize=24, fontweight="bold", pad=16)
    ax.tick_params(axis="both", labelsize=16, pad=2)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- assembly
def letter_font(height_px: int) -> ImageFont.FreeTypeFont:
    path = font_manager.findfont(font_manager.FontProperties(family=FONT, weight="bold"))
    size = 40
    while True:  # largest size whose "(b)" glyph box is <= the reference letter height
        f = ImageFont.truetype(path, size + 1)
        x0, y0, x1, y1 = f.getbbox("(b)")
        if y1 - y0 > height_px:
            return ImageFont.truetype(path, size)
        size += 1


def paste(canvas: Image.Image, im: Image.Image, box) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = box
    bw, bh = x1 - x0, y1 - y0
    s = min(bw / im.width, bh / im.height)
    w, h = round(im.width * s), round(im.height * s)
    im = im.resize((w, h), Image.LANCZOS)
    ox, oy = x0 + (bw - w) // 2, y0 + (bh - h) // 2
    canvas.paste(im, (ox, oy))
    return ox, oy, ox + w, oy + h


def vector_pdf(panels: dict, placed: dict, path: Path) -> None:
    """The composite as a vector PDF: each panel's PDF page, cropped to its ink, scaled into its placed box (the
    same boxes as the PNG) on a page of the canvas size at DPI_OUT; letters as TrueType text."""
    sys.path.insert(0, str(fig5_paths.FIG_DIR / ".pylib"))  # pypdf (pure Python), installed with pip --target
    from pypdf import PdfReader, PdfWriter, Transformation
    from matplotlib.textpath import TextPath

    pt = 72 / DPI_OUT
    page_w, page_h = CANVAS[0] * pt, CANVAS[1] * pt
    fig = plt.figure(figsize=(page_w / 72, page_h / 72), dpi=72)
    prop = font_manager.FontProperties(family=FONT, weight="bold")
    ext = TextPath((0, 0), "(b)", size=100, prop=prop).get_extents()
    size = LETTER_HEIGHT * pt / (ext.height / 100)
    for text, (x, y) in LETTERS.items():
        e = TextPath((0, 0), text, size=size, prop=prop).get_extents()
        fig.text(x * pt - e.x0, page_h - y * pt - e.y1, text, fontproperties=prop, fontsize=size, ha="left",
                 va="baseline", transform=None)
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf")
    plt.close(fig)
    writer = PdfWriter()
    base = writer.add_page(PdfReader(buf).pages[0])
    # 02b's expression-curve PDFs have a white background: merge them first so their padding sits under the rest
    for key, p in sorted(panels.items(), key=lambda kv: not kv[0].startswith("e_expr")):
        page = PdfReader(p["pdf"]).pages[0]
        mb = page.mediabox
        pw, ph, left, bottom = float(mb.width), float(mb.height), float(mb.left), float(mb.bottom)
        fx0, fy0, fx1, fy1 = p["frac"]
        cx0, cy0 = left + fx0 * pw, bottom + (1 - fy1) * ph  # ink crop, PDF coordinates (origin bottom left)
        cx1 = left + fx1 * pw
        ox0, _oy0, ox1, oy1 = placed[key]
        s = (ox1 - ox0) * pt / (cx1 - cx0)
        base.merge_transformed_page(page, Transformation().scale(s).translate(ox0 * pt - cx0 * s, page_h - oy1 * pt - cy0 * s))
    writer.add_metadata({"/Title": "Fig5 composite (Fig5_GRN, v4)"})
    with open(path, "wb") as fh:
        writer.write(fh)


def main() -> None:
    PANELS.mkdir(parents=True, exist_ok=True)
    style()
    strip.ordered_tfs_for_benchmark = v2.ordered_tfs_for_benchmark
    metrics = v2.filter_rows(pd.read_csv(v2.METRICS_TABLE))
    curves = filter_curve_rows(pd.read_csv(macro.CURVES_TABLE))
    pr_curves = pd.read_csv(fig5_paths.PUBLICATION_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_pr_curves.csv")

    print("== panels")
    panels = {
        "b_trrust": render_to_aspect(lambda fs: strip_panel(metrics, "TRRUST", fs), (10.6, 4.65), box_aspect("b_trrust"), "b_trrust"),
        "b_chip": render_to_aspect(lambda fs: strip_panel(metrics, "K562_ChIP", fs), (10.6, 4.65), box_aspect("b_chip"), "b_chip"),
        "c_trrust": render_to_aspect(lambda fs: macro_panel(curves, "TRRUST", fs), (4.0, 4.0), box_aspect("c_trrust"), "c_trrust"),
        "c_chip": render_to_aspect(lambda fs: macro_panel(curves, "K562_ChIP", fs), (4.0, 3.75), box_aspect("c_chip"), "c_chip"),
        "d": render_to_aspect(nes_panel, (6.6, 7.4), box_aspect("d"), "d"),
        "e_heat": render_to_aspect(heat_panel, (7.2, 3.55), box_aspect("e_heat"), "e_heat"),
    }
    for tf in ["RELA", "NFKB1", "JUN"]:
        panels[f"e_pr_{tf}"] = render_to_aspect(lambda fs, tf=tf: pr_panel(pr_curves, tf, fs), (4.3, 4.7),
                                                box_aspect(f"e_pr_{tf}"), f"e_pr_{tf}")
        stem = FIG5D / f"expr-curve_{tf}_chronopack_1lineage"  # 02b's PNG and PDF (same bounding box)
        panels[f"e_expr_{tf}"] = as_panel(Image.open(stem.with_suffix(".png")).convert("RGB"), stem.with_suffix(".pdf"))

    print("== composite")
    canvas = Image.new("RGB", CANVAS, "white")
    placed = {k: paste(canvas, p["im"], REF_BOXES[k]) for k, p in panels.items()}
    draw = ImageDraw.Draw(canvas)
    font = letter_font(LETTER_HEIGHT)
    for text, (x, y) in LETTERS.items():
        bx0, by0, _, _ = font.getbbox(text)
        draw.text((x - bx0, y - by0), text, font=font, fill="black")
    canvas.save(OUT / "Fig5_composite_v4.png", dpi=(DPI_OUT, DPI_OUT))
    vector_pdf(panels, placed, OUT / "Fig5_composite_v4.pdf")
    ref = Image.open(REF).convert("RGB").resize(CANVAS)
    Image.blend(ref, canvas, 0.5).save(OUT / "overlay_with_reference.png")  # layout check
    rows = [{"panel": k, "ref_box": REF_BOXES[k], "placed_box": v} for k, v in placed.items()]
    pd.DataFrame(rows).to_csv(OUT / "panel_boxes.csv", index=False)
    for r in rows:
        print(f"  {r['panel']:13s} ref {r['ref_box']} placed {r['placed_box']}")
    print(OUT / "Fig5_composite_v4.png")
    print("COMPOSITE_OK")


if __name__ == "__main__":
    main()
