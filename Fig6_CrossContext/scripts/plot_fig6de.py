#!/usr/bin/env python3
"""Fig6d/e (originally Fig7d/e): elastic-net regression screen of dynamic-mode kinetic outcomes.

6d  heatmap of the cross-validated R^2 lost when each feature family is dropped
    from the full model (cv_summary_results.csv, level ``drop_one_delta``;
    positive = the family carries unique predictive signal).
6e  dot plot of per-feature OLS partial R^2 from the adjusted individual-feature
    regressions (canonical9_ols_results.csv); dot area = partial R^2, colour = sign
    of the coefficient, black ring = FDR < 0.05. The ``cell_type`` adjustment
    covariate is not plotted.

Fig6_CrossContext_v2: both tables are the v4 kinetics tables of Fig6_Reproduction/reference/kinetics/
(which supersede v1's Run3 tables). They use the five disjoint families of the v4 regression screen,
and 6d's row labels are the table's ``family_label`` column.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "output" / ".mplconfig"))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fig6_common as combined  # noqa: E402

CV_INPUT = combined.INPUT_DIR / "cv_summary_results.csv"
OLS_INPUT = combined.INPUT_DIR / "canonical9_ols_results.csv"
OUT_D = combined.OUTPUT_DIR / "fig6d_drop_one_family_delta_r2"
OUT_E = combined.OUTPUT_DIR / "fig6e_individual_partial_r2"

# Kinetic outcomes, in column order for both panels.
OUTCOMES = [
    "baseline_exp",
    "t_onset",
    "t_peak",
    "peak_log2fc",
    "max_induction_rate",
    "owp",
    "owp_ddt",
    "endpoint_retention",
    "postpeak_drop",
    "max_decay_rate",
]
OUTCOME_LABELS = {
    "baseline_exp": "Baseline expr.",
    "t_onset": r"$t_{\mathrm{onset}}$",
    "t_peak": r"$t_{\mathrm{peak}}$",
    "peak_log2fc": "peak log2FC",
    "max_induction_rate": "max induction rate",
    "owp": "response persistence",
    "owp_ddt": "rate persistence",
    "endpoint_retention": "endpoint retention",
    "postpeak_drop": "post-peak drop",
    "max_decay_rate": "max decay rate",
}
# 6d rows (families, in order) and display labels; load_drop_one() checks the labels against the table's
# family_label column.
FAMILIES = [
    ("promoter_architecture", "Promoter architecture"),
    ("enhancer_architecture", "Enhancer architecture"),
    ("enhancer_burden", "Enhancer burden"),
    ("rna_stability", "RNA stability"),
    ("gene_length", "Gene length"),
]
# Two-line family labels for the narrow right margin of 6e in the composite.
FAMILY_WRAPPED = {
    "promoter_architecture": "Promoter\narchitecture",
    "enhancer_burden": "Enhancer\nburden",
    "rna_stability": "RNA\nstability",
    "gene_length": "Gene length",
}
# Individual features in row order, grouped by the family they belong to in 6d.
FEATURES = [
    ("RELA_in_PLS", "promoter_architecture"),
    ("AP-1_HOMER_in_PLS", "promoter_architecture"),
    ("CRE_HOMER_in_PLS", "promoter_architecture"),
    ("is_TATA_strict", "promoter_architecture"),
    ("median_log10_decay_rate", "rna_stability"),
    ("is_ARE", "rna_stability"),
    ("gene_length_kb", "gene_length"),
    ("n_pELS", "enhancer_burden"),
    ("n_dELS", "enhancer_burden"),
]
FEATURE_LABELS = {
    "RELA_in_PLS": "RELA",
    "AP-1_HOMER_in_PLS": "AP-1 (TRE)",
    "CRE_HOMER_in_PLS": "CRE",
    "is_TATA_strict": "TATA box",
    "median_log10_decay_rate": "decay rate",
    "is_ARE": "5'UTR ARE",
    "gene_length_kb": "gene length (kb)",
    "n_pELS": "# proximal enhancers",
    "n_dELS": "# distal enhancers",
}
FDR_CUTOFF = 0.05

CMAP = mpl.colormaps["RdBu_r"]
POS_COLOR = mpl.colors.to_hex(CMAP(0.86))
NEG_COLOR = mpl.colors.to_hex(CMAP(0.14))
INK = "#262626"
MUTED = "#6B6B6B"
BAND = "#F3F4F6"

SIZE_LEGEND = [0.01, 0.03, 0.07]

# Text and mark sizes for the standalone panels and the composite. 6e dot area (pt^2) is
# dot_floor + dot_slope * partial R^2, so area stays proportional above the floor.
STANDALONE = dict(
    tick=9.5, annot=7.6, cbar_tick=7.6, cbar_label=8.3, cell_gap=2.0, family=8,
    legend=8, legend_title=8.5, dot_floor=16.0, dot_slope=5100.0, ring=1.3, key_marker=7.5,
)
COMPOSITE = dict(
    tick=7, annot=5.6, cbar_tick=6, cbar_label=6.5, cell_gap=1.5, family=6.3,
    legend=6.3, legend_title=6.8, dot_floor=7.0, dot_slope=2700.0, ring=0.9, key_marker=5.5,
)
# Composite page size: 260 x 90 mm.
COMPOSITE_SIZE = (260 / 25.4, 90 / 25.4)
OUT_DE = combined.OUTPUT_DIR / "fig6de_composite"
# 6e row pitch relative to the original layout (0.8 = rows 20% closer together).
E_ROW_SCALE = 0.8


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 9,
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def fmt2(v: float) -> str:
    s = f"{v:.2f}"
    return "0.00" if s == "-0.00" else s


def check_family_labels(cv: pd.DataFrame) -> None:
    """The figure labels must be the table's family_label for every family shown."""
    table = cv.dropna(subset=["family_label"]).groupby("family")["family_label"].unique()
    bad = {f: (label, list(table.get(f, []))) for f, label in FAMILIES if list(table.get(f, [])) != [label]}
    if bad:
        raise ValueError(f"family labels differ from cv_summary_results.csv family_label: {bad}")


def load_level(level: str) -> pd.DataFrame:
    """Family x outcome matrix of r2_mean for one level of the CV summary (drop_one_delta, family_only)."""
    cv = pd.read_csv(CV_INPUT)
    check_family_labels(cv)
    cv = cv[cv["level"].eq(level)]
    mat = cv.pivot(index="family", columns="outcome", values="r2_mean")
    mat = mat.reindex(index=[f for f, _ in FAMILIES], columns=OUTCOMES)
    if mat.isna().any().any():
        raise ValueError(f"Missing {level} cells in CV summary")
    return mat


def load_drop_one() -> pd.DataFrame:
    return load_level("drop_one_delta")


def load_partial_r2() -> pd.DataFrame:
    ols = pd.read_csv(OLS_INPUT)
    features = pd.DataFrame(FEATURES, columns=["predictor", "family"])
    ols = ols.merge(features, on="predictor", how="inner")
    ols = ols[ols["outcome"].isin(OUTCOMES)].copy()
    if len(ols) != len(FEATURES) * len(OUTCOMES):
        raise ValueError(f"Expected one row per feature/outcome, got {len(ols)}")
    ols["x"] = ols["outcome"].map({o: i for i, o in enumerate(OUTCOMES)})
    ols["y"] = ols["predictor"].map({p: i for i, (p, _) in enumerate(FEATURES)})
    ols["direction"] = np.where(ols["coef"] >= 0, "positive", "negative")
    ols["significant"] = ols["padj_fdr"] < FDR_CUTOFF
    return ols.sort_values(["y", "x"])


def style_outcome_axis(ax: plt.Axes, sz: dict) -> None:
    ax.set_xticks(range(len(OUTCOMES)))
    ax.set_xticklabels([OUTCOME_LABELS[o] for o in OUTCOMES], rotation=35, ha="right", rotation_mode="anchor", fontsize=sz["tick"])
    ax.tick_params(axis="both", length=0, pad=0.4 * sz["tick"], labelsize=sz["tick"])
    for spine in ax.spines.values():
        spine.set_visible(False)


def save(fig: plt.Figure, out_prefix: Path) -> None:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".pdf"))
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300)
    plt.close(fig)
    print(f"Wrote {out_prefix.with_suffix('.pdf')}")
    print(f"Wrote {out_prefix.with_suffix('.png')}")


def write_values(df: pd.DataFrame, out_prefix: Path) -> None:
    values_csv = out_prefix.parent / f"{out_prefix.name}_values.csv"
    df.to_csv(values_csv, index=False)
    print(f"Wrote {values_csv}")


D_LABEL = "CV R² lost when family is dropped (ΔR²)"
D_LABEL_WRAPPED = "CV R² lost when\nfamily is dropped (ΔR²)"


def draw_d(
    ax: plt.Axes, cax: plt.Axes, mat: pd.DataFrame, sz: dict, orientation: str,
    row_labels: list[str] | None = None, label: str = D_LABEL, label_wrapped: str = D_LABEL_WRAPPED,
) -> None:
    """Annotated diverging heatmap (rows = mat's index, columns = OUTCOMES); also used by the family-only figure."""
    row_labels = [lab for _, lab in FAMILIES] if row_labels is None else row_labels
    n_rows = len(row_labels)
    vals = mat.to_numpy()
    vmax = np.ceil(np.abs(vals).max() / 0.02) * 0.02
    norm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
    ax.pcolormesh(
        np.arange(len(OUTCOMES) + 1) - 0.5,
        np.arange(n_rows + 1) - 0.5,
        vals,
        cmap=CMAP,
        norm=norm,
        edgecolors="white",
        linewidth=sz["cell_gap"],
    )
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v = vals[i, j]
            ax.text(
                j, i, fmt2(v),
                ha="center", va="center", fontsize=sz["annot"],
                color="white" if abs(v) > 0.55 * vmax else INK,
            )
    ax.set_xlim(-0.5, len(OUTCOMES) - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels)
    style_outcome_axis(ax, sz)

    cbar = ax.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation=orientation)
    cbar.set_ticks(np.linspace(-vmax, vmax, 5))
    cbar.ax.tick_params(labelsize=sz["cbar_tick"], length=2, color=MUTED, pad=2)
    cbar.outline.set_visible(False)
    if orientation == "vertical":
        cbar.set_label(label_wrapped, fontsize=sz["cbar_label"], labelpad=6)
    else:
        cax.xaxis.set_label_position("top")
        cbar.set_label(label, fontsize=sz["cbar_label"], labelpad=3)


def draw_e(
    ax: plt.Axes, pts: pd.DataFrame, sz: dict, legend_x: float, wrap_families: bool, key_y: float = 0.52,
) -> None:
    n_rows, n_cols = len(FEATURES), len(OUTCOMES)
    family_label = dict(FAMILIES)
    area = sz["dot_floor"] + sz["dot_slope"] * pts["partial_r2"]

    # Shade alternate feature families and label them in the right margin.
    groups = []
    for y, (_, fam) in enumerate(FEATURES):
        if groups and groups[-1][0] == fam:
            groups[-1][2] = y
        else:
            groups.append([fam, y, y])
    for k, (fam, y0, y1) in enumerate(groups):
        if k % 2 == 0:
            ax.axhspan(y0 - 0.5, y1 + 0.5, color=BAND, zorder=0, lw=0)
        label = FAMILY_WRAPPED[fam] if wrap_families else family_label[fam]
        ax.text(
            n_cols - 0.35, (y0 + y1) / 2, label,
            ha="left", va="center", fontsize=sz["family"], color=MUTED, style="italic",
            linespacing=1.1,
        )
    for x in range(n_cols):
        ax.axvline(x, color="#E2E2E2", lw=0.5, zorder=1)

    for sig in (False, True):
        keep = pts["significant"] == sig
        sub = pts[keep]
        ax.scatter(
            sub["x"], sub["y"], s=area[keep],
            c=np.where(sub["direction"] == "positive", POS_COLOR, NEG_COLOR),
            edgecolors="black" if sig else "white",
            linewidths=sz["ring"] if sig else 0.4,
            zorder=3 if sig else 2,
        )

    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([FEATURE_LABELS[p] for p, _ in FEATURES])
    style_outcome_axis(ax, sz)

    size_handles = [
        ax.scatter([], [], s=sz["dot_floor"] + sz["dot_slope"] * v, color="#9A9A9A", edgecolors="white", linewidths=0.4)
        for v in SIZE_LEGEND
    ]
    leg1 = ax.legend(
        size_handles, [f"{v:.2f}" for v in SIZE_LEGEND], title="Partial R²",
        loc="upper left", bbox_to_anchor=(legend_x, 1.0), frameon=False,
        fontsize=sz["legend"], title_fontsize=sz["legend_title"], labelspacing=1.1,
        borderpad=0, handletextpad=1.0, alignment="left",
    )
    ax.add_artist(leg1)
    m = sz["key_marker"]
    key_handles = [
        Line2D([0], [0], marker="o", ls="none", markersize=m, markerfacecolor=POS_COLOR,
               markeredgecolor="white", markeredgewidth=0.4, label="Positive"),
        Line2D([0], [0], marker="o", ls="none", markersize=m, markerfacecolor=NEG_COLOR,
               markeredgecolor="white", markeredgewidth=0.4, label="Negative"),
        Line2D([0], [0], marker="o", ls="none", markersize=m, markerfacecolor="white",
               markeredgecolor="black", markeredgewidth=sz["ring"], label=f"FDR < {FDR_CUTOFF:g}"),
    ]
    ax.legend(
        handles=key_handles, title="Coefficient", loc="upper left",
        bbox_to_anchor=(legend_x, key_y), frameon=False, fontsize=sz["legend"],
        title_fontsize=sz["legend_title"], borderpad=0, handletextpad=0.8, alignment="left",
    )


E_COLUMNS = ["outcome", "predictor", "family", "coef", "partial_r2", "pvalue", "padj_fdr", "significant"]


def plot_d(out_prefix: Path = OUT_D) -> None:
    mat = load_drop_one()
    fig, ax = plt.subplots(figsize=(8.4, 3.5))
    fig.subplots_adjust(left=0.245, right=0.85, top=0.96, bottom=0.27)
    cax = fig.add_axes([0.867, 0.35, 0.013, 0.58])
    draw_d(ax, cax, mat, STANDALONE, "vertical")
    save(fig, out_prefix)
    write_values(mat.stack().rename("delta_r2").reset_index(), out_prefix)


def plot_e(out_prefix: Path = OUT_E) -> None:
    pts = load_partial_r2()
    bottom_in, top_in = 1.10, 0.13
    axes_h = 3.168 * E_ROW_SCALE
    fig_h = bottom_in + axes_h + top_in
    fig, ax = plt.subplots(figsize=(9.6, fig_h))
    fig.subplots_adjust(left=0.20, right=0.70, top=1 - top_in / fig_h, bottom=bottom_in / fig_h)
    draw_e(ax, pts, STANDALONE, legend_x=1.22, wrap_families=True, key_y=0.40)
    save(fig, out_prefix)
    write_values(pts[E_COLUMNS], out_prefix)
    print(f"Points: {len(pts)}; FDR < {FDR_CUTOFF:g}: {int(pts['significant'].sum())}")


def plot_composite(out_prefix: Path = OUT_DE) -> None:
    """6d (left) and 6e (right) on one 260 x 90 mm landscape page, sharing an x-label baseline."""
    fig_w, fig_h = COMPOSITE_SIZE
    fig = plt.figure(figsize=COMPOSITE_SIZE)

    def inch_axes(x: float, y: float, w: float, h: float) -> plt.Axes:
        return fig.add_axes([x / fig_w, y / fig_h, w / fig_w, h / fig_h])

    bottom, top_pad = 0.74, 0.10
    e_h = (fig_h - bottom - top_pad) * E_ROW_SCALE
    d_cell_w, d_row_h = 0.285, 0.42
    e_cell_w = 0.26
    d_x, e_x = 1.42, 5.80
    d_h = d_row_h * len(FAMILIES)

    ax_d = inch_axes(d_x, bottom, d_cell_w * len(OUTCOMES), d_h)
    cax = inch_axes(d_x + 0.30 * d_cell_w * len(OUTCOMES), bottom + d_h + 0.30, 0.40 * d_cell_w * len(OUTCOMES), 0.07)
    draw_d(ax_d, cax, load_drop_one(), COMPOSITE, "horizontal")

    e_w = e_cell_w * len(OUTCOMES)
    ax_e = inch_axes(e_x, bottom, e_w, e_h)
    legend_x = (fig_w - 0.98 - e_x) / e_w
    draw_e(ax_e, load_partial_r2(), COMPOSITE, legend_x=legend_x, wrap_families=True, key_y=0.40)

    for letter, x in (("d", 0.06), ("e", e_x - 1.45)):
        fig.text(x / fig_w, 1 - 0.06 / fig_h, letter, ha="left", va="top", fontsize=10, fontweight="bold")
    save(fig, out_prefix)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel", choices=["d", "e", "de"], nargs="+", default=["d", "e", "de"],
                        help="d, e = standalone panels; de = 260 x 90 mm composite (default: all)")
    args = parser.parse_args()
    set_style()
    if "d" in args.panel:
        plot_d()
    if "e" in args.panel:
        plot_e()
    if "de" in args.panel:
        plot_composite()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
