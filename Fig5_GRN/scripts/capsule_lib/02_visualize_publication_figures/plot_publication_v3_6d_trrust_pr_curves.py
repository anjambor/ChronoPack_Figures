# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE = BASE_DIR / "benchmark_outputs" / "tables" / "publication_figures" / "system_broad_filtered_publication_v2_per_tf_pr_curves.csv"
OUT_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "5d"
SCENARIO = "system_broad_filtered_publication_v2"
METHOD_ORDER = [
    "ChronoPack_K562",
    "GENIE3",
    "Replogle2022_K562",
    "ChronoPack_THP1",
    "Yao2023_THP1",
]
METHOD_ORDER = fig5_paths.insert_prior(METHOD_ORDER, after="GENIE3")  # Fig5_GRN_v3 (with_prior)
METHOD_COLORS = {
    "ChronoPack_K562": "#D99000",
    "GENIE3": "#3D8B37",
    "Replogle2022_K562": "#CC79A7",
    "ChronoPack_THP1": "#D55E00",
    "Yao2023_THP1": "#3366AA",
    **fig5_paths.PRIOR_COLORS,
}
METHOD_LABELS = {
    "ChronoPack_K562": "Chrono (K562)",
    "GENIE3": "GENIE3 (K562)",
    "Replogle2022_K562": "Replogle (K562)",
    "ChronoPack_THP1": "Chrono (THP1)",
    "Yao2023_THP1": "Yao (THP1)",
    **{m: f"{label} (K562)" for m, label in fig5_paths.PRIOR_LABELS.items()},
}
YLIMS = {
    "RELA": 0.7,
    "NFKB1": 0.7,
    "JUN": 0.10,
}
VARIANT_YLIM = 1.0
TITLE_SIZE = 16
LABEL_SIZE = 15
TICK_SIZE = 15
LEGEND_SIZE = 9


def save(fig: plt.Figure, stem: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_tf(curves: pd.DataFrame, tf: str, *, ymax: float, suffix: str = "") -> None:
    panel = curves[
        curves["scenario"].eq(SCENARIO)
        & curves["benchmark"].eq("TRRUST")
        & curves["scope"].eq("per_tf")
        & curves["TF"].eq(tf)
        & curves["method"].isin(METHOD_ORDER)
    ].copy()

    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    for method in METHOD_ORDER:
        d = panel[panel["method"].eq(method)].sort_values("k")
        if d.empty:
            continue
        observed = d[~d["is_expectation"]].dropna(subset=["recall", "precision"])
        expected = d[d["is_expectation"]].dropna(subset=["recall", "precision"])
        color = METHOD_COLORS[method]
        if len(observed):
            ax.plot(
                observed["recall"],
                observed["precision"],
                drawstyle="steps-post",
                lw=2.2,
                color=color,
                label=METHOD_LABELS[method],
            )
        if len(expected):
            start = observed.tail(1) if len(observed) else expected.head(1)
            tail = pd.concat([start, expected]).drop_duplicates(["recall", "precision"])
            ax.plot(tail["recall"], tail["precision"], lw=1.9, ls="--", color=color, alpha=0.75)

    ax.set_xlim(0, 0.25)
    ax.set_ylim(0, ymax)
    if tf == "JUN" and suffix == "":
        ax.set_yticks(np.arange(0, ymax + 0.001, 0.04))
    ax.set_box_aspect(1)
    ax.set_xlabel("Recall", fontsize=LABEL_SIZE)
    ax.set_ylabel("Precision", fontsize=LABEL_SIZE)
    ax.set_title(tf, fontsize=TITLE_SIZE, pad=10)
    ax.tick_params(axis="both", labelsize=TICK_SIZE)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if tf == "RELA":
        ax.legend(frameon=False, fontsize=LEGEND_SIZE, loc="upper right", handlelength=1.6, borderaxespad=0.2)
    fig.tight_layout()
    stem = f"publication_pr_TRRUST_{tf}{suffix}"
    save(fig, stem)
    print(OUT_DIR / f"{stem}.pdf")
    print(OUT_DIR / f"{stem}.png")


def main() -> None:
    curves = pd.read_csv(TABLE)
    for tf in ["RELA", "NFKB1", "JUN"]:
        plot_tf(curves, tf, ymax=YLIMS[tf])
        plot_tf(curves, tf, ymax=VARIANT_YLIM, suffix="_ymax1p0")


if __name__ == "__main__":
    main()
