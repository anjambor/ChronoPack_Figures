# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from generate_system_broad_filtered_assets import (
    CMAP,
    SYSTEM_BLUE,
    apply_annotation_colors,
    format_annotations,
    normalize_over_method_rows,
    style_heatmap_axis,
)
from generate_system_broad_filtered_publication_figures import (
    PR_METHOD_LABELS,
    PR_METHOD_ORDER,
    add_compact_row_label_bands,
    add_trrust_tf_cell_type_row_groups,
)


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables" / "publication_figures"
OUT_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "5d"

METRICS_TABLE = TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv"
SOURCE_OUT = TABLE_DIR / "system_broad_filtered_publication_v3_trrust_core_tf_metrics_reference_palette.csv"
CORE_TFS = ["RELA", "NFKB1", "JUN"]
OUT_STEM = "publication_heatmap_TRRUST_RELA_NFKB1_JUN_EPR_AUPRCr"
FIGSIZE = (8.0, 4.0)


def trrust_core_tf_matrix(metrics: pd.DataFrame) -> pd.DataFrame:
    sub = metrics[
        metrics["benchmark"].eq("TRRUST")
        & metrics["scope"].eq("per_tf")
        & metrics["method"].isin(PR_METHOD_ORDER)
        & metrics["TF"].isin(CORE_TFS)
    ]
    cols = []
    for tf in CORE_TFS:
        tf_sub = sub[sub["TF"].eq(tf)].set_index("method")
        cols.append(tf_sub["EPR"].reindex(PR_METHOD_ORDER).rename((tf, "EPR")))
        cols.append(tf_sub["AUPRC_ratio"].reindex(PR_METHOD_ORDER).rename((tf, "AUPRCr")))
    mat = pd.concat(cols, axis=1).rename(index=PR_METHOD_LABELS)
    mat.columns = pd.MultiIndex.from_tuples(mat.columns)
    return mat


def save(fig: plt.Figure) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{OUT_STEM}.png", dpi=300, bbox_inches=None)
    fig.savefig(OUT_DIR / f"{OUT_STEM}.pdf", dpi=300, bbox_inches=None)
    plt.close(fig)


def main() -> None:
    sns.set_theme(style="white", font_scale=0.9)
    mat = trrust_core_tf_matrix(pd.read_csv(METRICS_TABLE))
    mat.to_csv(SOURCE_OUT)
    color = normalize_over_method_rows(mat)

    rows = len(mat.index)  # Fig5_GRN_v3: v2 has 5 rows; with_prior adds 3
    fig, ax = plt.subplots(figsize=FIGSIZE if rows <= 5 else (FIGSIZE[0], FIGSIZE[1] * rows / 5))
    sns.heatmap(
        color,
        annot=format_annotations(mat),
        fmt="",
        cmap=CMAP,
        vmin=0,
        vmax=1,
        cbar=False,
        linewidths=1.25,
        linecolor="white",
        ax=ax,
        annot_kws={"fontsize": 10.5},
    )
    apply_annotation_colors(ax, color)
    add_compact_row_label_bands(ax, [str(x) for x in mat.index], width=0.18, fontsize=10.5)
    add_trrust_tf_cell_type_row_groups(ax)
    style_heatmap_axis(ax, "")
    ax.set_xticklabels(mat.columns.get_level_values(1), rotation=0, fontsize=10.5)
    for tf in CORE_TFS:
        idxs = [i for i, val in enumerate(mat.columns.get_level_values(0)) if val == tf]
        center = (min(idxs) + max(idxs) + 1) / 2
        ax.text(
            center / mat.shape[1],
            1.20,
            tf,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=12,
            fontweight="bold",
            color=SYSTEM_BLUE,
        )
    for x in range(2, mat.shape[1], 2):
        ax.axvline(x, color="white", lw=5.0)
    scale = min(1.0, 5 / rows)  # Fig5_GRN_v3: keep v2's absolute top / bottom margins on a taller figure
    fig.subplots_adjust(left=0.28, right=0.98, top=1 - 0.38 * scale, bottom=0.08 * scale)
    save(fig)
    print(OUT_DIR / f"{OUT_STEM}.pdf")
    print(OUT_DIR / f"{OUT_STEM}.png")
    print(SOURCE_OUT)


if __name__ == "__main__":
    main()
