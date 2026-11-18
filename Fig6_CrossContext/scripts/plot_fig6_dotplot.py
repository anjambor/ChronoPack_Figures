#!/usr/bin/env python3
"""Fig6 (originally Fig7): transposed combined ORA and PLS motif/covariate dynamic-mode dot plot.

Adapted from
/new-stg/home/alex/tmp/fig3_dm_dotplot/ora_v3/analysis/plot_combined_ora_pls_motif_dotplot_transposed.py.
Writes one figure per covariate (logistic regression screen) variant.
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

# Only the v2 (percentile-covariate) screen is copied into input/.
VARIANTS = {
    "v2": (
        combined.INPUT_DIR / "logistic_regression_screen_results_v3_b3_pctl_fg4.csv",
        combined.OUTPUT_DIR / "fig6_ora_cov_pls_motifs_transposed_v2_pctl_fg4",
    ),
}
# Outline cutoffs (ORA q, motif q, covariate p_adj) and output-name suffix.
# "mixed" is the reference figure's convention.
CUTOFF_SETS = {
    "mixed": (combined.ORA_Q_CUTOFF, combined.MOTIF_Q_CUTOFF, combined.COVARIATE_PADJ_CUTOFF, ""),
    "q05": (0.05, 0.05, 0.05, "_q05"),
    "q10": (0.10, 0.10, 0.10, "_q10"),
}
POINT_COLUMNS = [
    "panel", "display", "cell", "cluster", "effect", "p_value",
    "plot_neg_log10_p", "q_value", "outlined",
]
ORA_AREA_SCALE = 0.80
MODE_EDGE_MARGIN = 0.15
CELL_SPACING = 0.24
MODE_PANEL_WIDTH = 2 * MODE_EDGE_MARGIN + CELL_SPACING
MODE_GAP = 0.08
CELL_OFFSETS = {
    "K562a": MODE_EDGE_MARGIN / MODE_PANEL_WIDTH,
    "THP1": (MODE_EDGE_MARGIN + CELL_SPACING) / MODE_PANEL_WIDTH,
}


def ora_dot_area(gene_ratio: float | None) -> float:
    if gene_ratio is None or not np.isfinite(gene_ratio):
        return 0.0
    return 16.0 + ORA_AREA_SCALE * 470.0 * (
        max(0.0, min(gene_ratio, 0.7)) / 0.7
    ) ** 0.5


def assign_context_x() -> tuple[dict[tuple[str, str], float], pd.DataFrame]:
    records = []
    lookup = {}
    for mode_idx, mode in enumerate(combined.MODE_LEVELS):
        mode_left = mode_idx * (MODE_PANEL_WIDTH + MODE_GAP)
        mode_right = mode_left + MODE_PANEL_WIDTH
        for cell in combined.CELL_LEVELS:
            x = mode_left + CELL_OFFSETS[cell] * MODE_PANEL_WIDTH
            lookup[(cell, mode)] = x
            records.append(
                {
                    "mode": mode,
                    "cell": cell,
                    "x": x,
                    "display": "K562" if cell == "K562a" else "THP1",
                    "mode_left": mode_left,
                    "mode_right": mode_right,
                    "mode_xmid": (mode_left + mode_right) / 2,
                }
            )
    return lookup, pd.DataFrame(records)


def assign_feature_y(ora_terms: pd.DataFrame, motif_rows: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    blocks = []
    y = 0.0
    row_step = 0.46
    internal_gap = 0.12
    separation_gap = 0.44

    ora_terms = ora_terms.copy()
    for group, sub in ora_terms.groupby("mode", sort=False):
        idx = sub.index
        ys = y + np.arange(len(sub)) * row_step
        ora_terms.loc[idx, "feature_y"] = ys
        blocks.append(
            {
                "frame": "ORA",
                "group": group,
                "ymin": ys.min() - row_step / 2,
                "ymax": ys.max() + row_step / 2,
                "ymid": ys.mean(),
            }
        )
        y = ys.max() + row_step + internal_gap

    y = blocks[-1]["ymax"] + separation_gap + row_step / 2
    motif_rows = motif_rows.copy()
    for group, sub in motif_rows.groupby("group", sort=False):
        idx = sub.index
        ys = y + np.arange(len(sub)) * row_step
        motif_rows.loc[idx, "feature_y"] = ys
        blocks.append(
            {
                "frame": group,
                "group": group,
                "ymin": ys.min() - row_step / 2,
                "ymax": ys.max() + row_step / 2,
                "ymid": ys.mean(),
            }
        )
        gap = separation_gap if group == "Covariates" else internal_gap
        y = ys.max() + row_step + gap

    return ora_terms, motif_rows, pd.DataFrame(blocks)


def draw_background(
    ax,
    feature_blocks: pd.DataFrame,
    context_cols: pd.DataFrame,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
) -> None:
    block_colors = ["#FAFAFA", "#F2F5F7"]
    for _, frame_blocks in feature_blocks.groupby("frame", sort=False):
        ax.add_patch(
            mpl.patches.Rectangle(
                (x_min, frame_blocks["ymin"].min()),
                x_max - x_min,
                frame_blocks["ymax"].max() - frame_blocks["ymin"].min(),
                facecolor="#FAFAFA",
                edgecolor="none",
                zorder=0,
            )
        )
    for i, block in feature_blocks.iterrows():
        color = "#FAFAFA" if block["group"] == "Covariates" else block_colors[i % 2]
        ax.add_patch(
            mpl.patches.Rectangle(
                (x_min, block["ymin"]),
                x_max - x_min,
                block["ymax"] - block["ymin"],
                facecolor=color,
                edgecolor="none",
                zorder=0,
            )
        )
    for edge in feature_blocks["ymax"].iloc[:-1]:
        ax.plot([x_min, x_max], [edge, edge], color="#BFBFBF", linewidth=0.8, zorder=1)

    for mode, sub in context_cols.groupby("mode", sort=False):
        xmin = sub["mode_left"].iloc[0]
        xmax = sub["mode_right"].iloc[0]
        ax.add_patch(
            mpl.patches.Rectangle(
                (xmin, y_min),
                xmax - xmin,
                y_max - y_min,
                facecolor="#F8F8F8",
                edgecolor="none",
                zorder=-1,
            )
        )
        ax.plot([xmin, xmin], [y_min, y_max], color="#D2D2D2", linewidth=0.6, zorder=1)
        ax.text(
            sub["mode_xmid"].iloc[0],
            y_min - 0.48,
            combined.MODE_LABELS[mode],
            ha="center",
            va="top",
            fontsize=9,
            fontweight="bold",
        )
        for _, row in sub.iterrows():
            ax.text(
                row["x"],
                y_min - 0.12,
                row["display"],
                ha="center",
                va="center",
                fontsize=9,
            )
    right_edge = context_cols["mode_right"].max()
    ax.plot([right_edge, right_edge], [y_min, y_max], color="#D2D2D2", linewidth=0.6, zorder=1)


def plot_variant(covariate_csv: Path, out_prefix: Path, cutoff_set: str = "mixed") -> None:
    ora_q_cutoff, motif_q_cutoff, cov_padj_cutoff, suffix = CUTOFF_SETS[cutoff_set]
    out_prefix = out_prefix.with_name(out_prefix.name + suffix)
    if ora_q_cutoff == motif_q_cutoff == cov_padj_cutoff:
        outline_label = f"outlined: q < {ora_q_cutoff:g}"
    else:
        outline_label = "outlined: q/padj cutoff"
    context_lookup, context_cols = assign_context_x()
    ora_terms = combined.load_ora_terms()
    motif_rows = combined.load_motif_rows()
    ora_terms, motif_rows, feature_blocks = assign_feature_y(ora_terms, motif_rows)

    y_lookup, _, _ = combined.make_row_geometry()
    ora_plot = combined.load_ora_plot(ora_terms, y_lookup).merge(
        ora_terms[["display", "feature_y"]],
        on="display",
        how="left",
    )
    motif_plot, cov_plot = combined.load_motif_plot(
        y_lookup,
        covariate_csv,
        motif_q_cutoff=motif_q_cutoff,
        covariate_padj_cutoff=cov_padj_cutoff,
    )
    motif_plot = motif_plot.merge(motif_rows[["display", "feature_y"]], on="display", how="left")
    cov_plot = cov_plot.merge(motif_rows[["display", "feature_y"]], on="display", how="left")

    for df in [ora_plot, motif_plot, cov_plot]:
        df["x"] = df.apply(lambda r: context_lookup[(r["cell"], r["cluster"])], axis=1)

    feature_labels = pd.concat(
        [
            ora_terms[["display", "feature_y"]],
            motif_rows[["display", "feature_y"]],
        ],
        ignore_index=True,
    )
    x_min = context_cols["mode_left"].min()
    x_max = context_cols["mode_right"].max()
    y_min = feature_blocks["ymin"].min()
    y_max = feature_blocks["ymax"].max()

    fig, ax = plt.subplots(figsize=(11.76, 12.8))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    draw_background(ax, feature_blocks, context_cols, x_min, x_max, y_min, y_max)

    for y in feature_labels["feature_y"]:
        ax.plot([x_min, x_max], [y, y], color="#ECECEC", linewidth=0.45, linestyle=":", zorder=1)

    ora_norm = mpl.colors.TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=5.0)
    ora_cmap = mpl.colormaps["RdBu_r"]
    motif_norm = mpl.colors.Normalize(vmin=-2.0, vmax=2.0, clip=True)
    motif_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "motif_log2_enrichment",
        ["#1b9e77", "#f7f7f7", "#d95f02"],
    )
    cov_norm = mpl.colors.Normalize(vmin=-2.0, vmax=2.0, clip=True)
    cov_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "covariate_log2_odds",
        ["#6b3f8f", "#f7f7f7", "#df8f2d"],
    )

    missing = ora_plot[~ora_plot["present"]]
    if not missing.empty:
        ax.scatter(missing["x"], missing["feature_y"], marker="+", s=20, color="#CFCFCF", linewidths=0.7, zorder=2)

    for cell, marker in {"K562a": "o", "THP1": "s"}.items():
        sub = ora_plot[ora_plot["present"] & ora_plot["cell"].eq(cell)]
        sig = sub[sub["qvalue"] < ora_q_cutoff]
        nonsig = sub[~(sub["qvalue"] < ora_q_cutoff)]
        if not nonsig.empty:
            ax.scatter(
                nonsig["x"],
                nonsig["feature_y"],
                s=nonsig["GeneRatio_val"].map(ora_dot_area),
                c=nonsig["log2FE"],
                cmap=ora_cmap,
                norm=ora_norm,
                marker=marker,
                edgecolors="none",
                zorder=3,
            )
        if not sig.empty:
            ax.scatter(
                sig["x"],
                sig["feature_y"],
                s=sig["GeneRatio_val"].map(ora_dot_area),
                c=sig["log2FE"],
                cmap=ora_cmap,
                norm=ora_norm,
                marker=marker,
                edgecolors="black",
                linewidths=0.62,
                zorder=4,
            )

    for _, row in cov_plot.iterrows():
        ax.scatter(
            row["x"],
            row["feature_y"],
            s=row["area"],
            c=[cov_cmap(cov_norm(row["effect"]))],
            marker="o" if row["cell"] == "K562a" else "s",
            edgecolors="black" if row["outlined"] else "none",
            linewidths=0.7 if row["outlined"] else 0,
            zorder=4,
        )
    for _, row in motif_plot.iterrows():
        ax.scatter(
            row["x"],
            row["feature_y"],
            s=row["area"],
            c=[motif_cmap(motif_norm(row["effect"]))],
            marker="o" if row["cell"] == "K562a" else "s",
            edgecolors="black" if row["outlined"] else "none",
            linewidths=0.7 if row["outlined"] else 0,
            zorder=4,
        )

    for _, frame_blocks in feature_blocks.groupby("frame", sort=False):
        ax.add_patch(
            mpl.patches.Rectangle(
                (x_min, frame_blocks["ymin"].min()),
                x_max - x_min,
                frame_blocks["ymax"].max() - frame_blocks["ymin"].min(),
                fill=False,
                linewidth=0.8,
                edgecolor="#4D4D4D",
                zorder=5,
            )
        )

    ax.set_xticks([])
    ax.set_yticks([])
    for _, row in feature_labels.iterrows():
        ax.text(
            x_min - 0.045,
            row["feature_y"],
            row["display"],
            ha="right",
            va="center",
            fontsize=9,
            clip_on=False,
        )
    ax.invert_yaxis()
    ax.set_xlim(x_min - 0.28, x_max + 0.08)
    ax.set_ylim(y_max + 0.35, y_min - 0.68)

    for spine in ax.spines.values():
        spine.set_visible(False)

    size_handles = [
        plt.scatter([], [], s=ora_dot_area(v), marker="o", facecolor="#8F8F8F", edgecolor="#8F8F8F")
        for v in [0.05, 0.10, 0.25, 0.50]
    ]
    leg1 = ax.legend(size_handles, ["0.05", "0.10", "0.25", "0.50"], title="ORA GeneRatio", loc="lower center", bbox_to_anchor=(0.22, -0.12), ncol=4, frameon=False, fontsize=8, title_fontsize=8.5)
    ax.add_artist(leg1)
    nlp_handles = [plt.scatter([], [], s=combined.motif_area(v), marker="o", facecolor="#8F8F8F", edgecolor="#8F8F8F") for v in [1, 2, 3]]
    leg2 = ax.legend(nlp_handles, ["1", "2", ">=3"], title="Motif/cov -log10(p)", loc="lower center", bbox_to_anchor=(0.63, -0.12), ncol=3, frameon=False, fontsize=8, title_fontsize=8.5)
    ax.add_artist(leg2)
    shape_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#8F8F8F", markeredgecolor="#8F8F8F", markersize=6.8, label="K562"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#8F8F8F", markeredgecolor="#8F8F8F", markersize=6.8, label="THP1"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="black", markeredgewidth=0.8, markersize=6.8, label=outline_label),
    ]
    leg3 = ax.legend(handles=shape_handles, loc="lower center", bbox_to_anchor=(0.50, -0.165), ncol=3, frameon=False, fontsize=8, handletextpad=0.45, columnspacing=0.85)
    ax.add_artist(leg3)

    cbar_ax1 = fig.add_axes([0.16, 0.055, 0.20, 0.012])
    fig.colorbar(mpl.cm.ScalarMappable(norm=ora_norm, cmap=ora_cmap), cax=cbar_ax1, orientation="horizontal").set_label("ORA log2(Fold enrichment)", fontsize=8.3, labelpad=2)
    cbar_ax2 = fig.add_axes([0.43, 0.055, 0.18, 0.012])
    fig.colorbar(mpl.cm.ScalarMappable(norm=motif_norm, cmap=motif_cmap), cax=cbar_ax2, orientation="horizontal").set_label("PLS motif log2 enrichment", fontsize=8.3, labelpad=2)
    cbar_ax3 = fig.add_axes([0.68, 0.055, 0.18, 0.012])
    fig.colorbar(mpl.cm.ScalarMappable(norm=cov_norm, cmap=cov_cmap), cax=cbar_ax3, orientation="horizontal").set_label("Covariate log2 odds", fontsize=8.3, labelpad=2)
    for cax in [cbar_ax1, cbar_ax2, cbar_ax3]:
        cax.tick_params(labelsize=7.4, length=2)

    fig.subplots_adjust(left=0.32, right=0.985, top=0.88, bottom=0.21)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".pdf"))
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300)
    plt.close(fig)

    ora_points = ora_plot[ora_plot["present"]].assign(
        panel="ORA",
        effect=ora_plot["log2FE"],
        p_value=np.nan,
        plot_neg_log10_p=np.nan,
        q_value=ora_plot["qvalue"],
        outlined=ora_plot["qvalue"] < ora_q_cutoff,
    )
    points = pd.concat(
        [
            ora_points[POINT_COLUMNS],
            cov_plot.assign(panel="Covariates")[POINT_COLUMNS],
            motif_plot.assign(panel="PLS motifs")[POINT_COLUMNS],
        ],
        ignore_index=True,
    )
    points_csv = out_prefix.parent / f"{out_prefix.name}_points.csv"
    points.to_csv(points_csv, index=False)

    print(f"Covariates: {covariate_csv.name}")
    print(f"Wrote {out_prefix.with_suffix('.pdf')}")
    print(f"Wrote {out_prefix.with_suffix('.png')}")
    print(f"Wrote {points_csv}")
    print(f"Feature rows: {len(feature_labels)}")
    print(f"Dynamic-mode columns: {len(context_cols)}")
    print("Outlined points per panel:", points.groupby("panel", sort=False)["outlined"].sum().to_dict())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--variant",
        choices=sorted(VARIANTS),
        nargs="+",
        default=["v2"],
        help="covariate variant(s) to plot (default: v2, the variant used for Fig6)",
    )
    parser.add_argument(
        "--cutoffs",
        choices=list(CUTOFF_SETS),
        nargs="+",
        default=["q05"],
        help="outline cutoff set(s): mixed = ORA q<0.05, motif/covariate q<0.10 (default: q05, used for Fig6)",
    )
    args = parser.parse_args()
    for name in args.variant:
        covariate_csv, out_prefix = VARIANTS[name]
        for cutoff_set in args.cutoffs:
            print(f"== {name} / {cutoff_set}")
            plot_variant(covariate_csv, out_prefix, cutoff_set)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
