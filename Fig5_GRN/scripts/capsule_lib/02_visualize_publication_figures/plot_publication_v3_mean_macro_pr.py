# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd

import benchmark_grn_systems as b
from generate_system_broad_filtered_assets import BENCHMARK_LABELS
from generate_system_broad_filtered_publication_figures import PR_METHOD_ORDER


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables" / "publication_figures"
OUT_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "5c"

CURVES_TABLE = TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_pr_curves.csv"
YLIMS = {
    "TRRUST": (0.0, 0.35),
    "K562_ChIP": (0.0, 0.5),
}
YTICKS = {  # Fig5 addition: default y ticks (02_reproduce_figures overrides YLIMS / YTICKS for FIG5_UNIVERSE=shared)
    "TRRUST": [0.0, 0.1, 0.2, 0.3],
    "K562_ChIP": [0.0, 0.2, 0.4],
}
TITLE_SIZE = 16
LABEL_SIZE = 15
TICK_SIZE = 15


def save(fig: plt.Figure, stem: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_segmented_line(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    dashed: np.ndarray,
    *,
    color: str,
    label: str,
) -> None:
    label_used = False
    for target_dashed, linestyle, linewidth, alpha in [
        (False, "-", 2.4, 1.0),
        (True, "--", 2.0, 0.72),
    ]:
        mask = dashed == target_dashed
        if not mask.any():
            continue
        idx = np.flatnonzero(mask)
        splits = np.where(np.diff(idx) > 1)[0] + 1
        for group in np.split(idx, splits):
            if len(group) == 1:
                group = np.array([max(0, group[0] - 1), group[0]])
            ax.plot(
                x[group],
                y[group],
                lw=linewidth,
                ls=linestyle,
                color=color,
                alpha=alpha,
                label=label if not label_used else None,
            )
            label_used = True


def plot_macro_pr_for_benchmark(
    curves: pd.DataFrame,
    benchmark: str,
    *,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    out_suffix: str = "",
    yticks: list[float] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    recall_grid = np.linspace(0, 1, 250)
    for method in PR_METHOD_ORDER:
        method_curves = []
        method_completion_masks = []
        panel = curves[(curves["benchmark"].eq(benchmark)) & (curves["method"].eq(method))]
        for tf in sorted(panel["TF"].dropna().unique()):
            d = panel[panel["TF"].eq(tf)].sort_values(["recall"])
            if d.empty:
                continue
            d = d.dropna(subset=["recall", "precision"]).drop_duplicates("recall", keep="last")
            x = d["recall"].to_numpy(float)
            y = d["precision"].to_numpy(float)
            if len(x) == 0:
                continue
            order = np.argsort(x)
            method_curves.append(np.interp(recall_grid, x[order], y[order], left=y[order][0], right=y[order][-1]))
            observed = d[~d["is_expectation"]]
            max_observed_recall = observed["recall"].max() if len(observed) else -np.inf
            method_completion_masks.append(recall_grid > max_observed_recall)
        if not method_curves:
            continue
        mean_precision = np.nanmean(np.vstack(method_curves), axis=0)
        completion_mask = np.any(np.vstack(method_completion_masks), axis=0)
        plot_segmented_line(
            ax,
            recall_grid,
            mean_precision,
            completion_mask,
            color=b.METHOD_COLORS[method],
            label=b.METHOD_LABELS[method],
        )

    ax.set_xlim(*(xlim if xlim is not None else (0, 0.25 if benchmark == "TRRUST" else 1)))
    ax.set_ylim(*(ylim if ylim is not None else YLIMS[benchmark]))
    if benchmark == "TRRUST":
        ax.set_yticks(YTICKS["TRRUST"])
    else:
        if xlim is not None and np.isclose(xlim[1], 0.5):
            ax.set_xticks([0.0, 0.5])
        else:
            ax.set_xticks([0.0, 0.5, 1.0])
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        if yticks is not None:  # Fig5 addition (K562 zoom variants)
            ax.set_yticks(yticks)
        elif ylim is not None and ylim[0] > 0:
            ax.set_yticks([tick for tick in [0.2, 0.4] if ylim[0] <= tick <= ylim[1]])
        else:
            ax.set_yticks(YTICKS["K562_ChIP"])
    ax.set_box_aspect(1)
    ax.set_xlabel("Recall", fontsize=LABEL_SIZE)
    ax.set_ylabel("Precision", fontsize=LABEL_SIZE)
    ax.set_title(BENCHMARK_LABELS[benchmark], fontsize=TITLE_SIZE, pad=10)
    ax.tick_params(axis="both", labelsize=TICK_SIZE)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if benchmark == "TRRUST":
        ax.legend(frameon=False, fontsize=6.5, loc="upper right", handlelength=1.6, borderaxespad=0.2)
    fig.tight_layout()
    save(fig, f"publication_mean_macro_pr_{benchmark}{out_suffix}")
    print(OUT_DIR / f"publication_mean_macro_pr_{benchmark}{out_suffix}.pdf")
    print(OUT_DIR / f"publication_mean_macro_pr_{benchmark}{out_suffix}.png")


def main() -> None:
    curves = pd.read_csv(CURVES_TABLE)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_macro_pr_for_benchmark(curves, benchmark)


if __name__ == "__main__":
    main()
