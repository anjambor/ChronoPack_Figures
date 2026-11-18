#!/usr/bin/env python3
"""Full-model cross-validated R^2 per kinetic outcome (publication version of refs/figureA_full_model_cv_r2.png).

Bars: elastic-net full-model CV R^2 (cv_summary_results.csv, level ``full_model``, family ``full``), mean +- SD
over the 10 CV repeats; bar colour = outcome group (as in the reference).
Diamonds: mean CV R^2 of the cell-type-only model (cv_celltype_only_results.csv; held-out intercept + cell-type
fits on the same outcome-specific folds). Both are Fig6_Reproduction's v4 kinetics tables.

Outcome order and labels are those of 6d/6e (plot_fig6de.py), except onset / peak time spelled out. Writes output/fig6_full_model_cv_r2.{pdf,png}
and a _values.csv with the plotted numbers (including the full-model min / max over CV repeats and the
cell-type-only SD, which are not drawn).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "output" / ".mplconfig"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fig6_common as combined  # noqa: E402
import plot_fig6de as de  # noqa: E402

CELLTYPE_INPUT = combined.INPUT_DIR / "cv_celltype_only_results.csv"
OUT = combined.OUTPUT_DIR / "fig6_full_model_cv_r2"

# Outcome groups and colours (grouping as in the reference figure).
GROUPS = {
    "Baseline": ["baseline_exp"],
    "Early": ["t_onset", "t_peak", "peak_log2fc", "max_induction_rate"],
    "Persistence": ["owp", "owp_ddt", "endpoint_retention", "postpeak_drop"],
    "Decline": ["max_decay_rate"],
}
GROUP_COLORS = {"Baseline": "#A3A3A3", "Early": "#4F86C0", "Persistence": "#7DA75A", "Decline": "#E39A45"}
OUTCOME_GROUP = {o: g for g, outcomes in GROUPS.items() for o in outcomes}

SIZE = (120 / 25.4, 72 / 25.4)  # 120 x 72 mm
FS = dict(tick=7.5, label=8, legend=7)
GRID = "#E6E6E6"
# 6d/6e outcome labels, with plain-word timing labels for this figure.
OUTCOME_LABELS = {**de.OUTCOME_LABELS, "t_onset": "onset time", "t_peak": "peak time"}


def load_values() -> pd.DataFrame:
    cv = pd.read_csv(de.CV_INPUT)
    full = cv[(cv["level"] == "full_model") & (cv["family"] == "full")].set_index("outcome")
    ct = pd.read_csv(CELLTYPE_INPUT).set_index("outcome")
    missing = [o for o in de.OUTCOMES if o not in full.index or o not in ct.index]
    if missing:
        raise ValueError(f"outcomes missing from the inputs: {missing}")
    out = pd.DataFrame({
        "outcome": de.OUTCOMES,
        "label": [OUTCOME_LABELS[o] for o in de.OUTCOMES],
        "group": [OUTCOME_GROUP[o] for o in de.OUTCOMES],
        "r2_full_mean": full.loc[de.OUTCOMES, "r2_mean"].to_numpy(),
        "r2_full_sd": full.loc[de.OUTCOMES, "r2_sd"].to_numpy(),
        "r2_full_min": full.loc[de.OUTCOMES, "r2_min"].to_numpy(),
        "r2_full_max": full.loc[de.OUTCOMES, "r2_max"].to_numpy(),
        "n_repeats": full.loc[de.OUTCOMES, "n_repeats"].to_numpy(),
        "r2_celltype_only": ct.loc[de.OUTCOMES, "r2_mean"].to_numpy(),
        "r2_celltype_only_sd": ct.loc[de.OUTCOMES, "r2_sd"].to_numpy(),
    })
    out["r2_gain_over_celltype"] = out["r2_full_mean"] - out["r2_celltype_only"]
    return out


def plot(values: pd.DataFrame, out_prefix: Path = OUT) -> None:
    fig, ax = plt.subplots(figsize=SIZE)
    fig.subplots_adjust(left=0.115, right=0.985, bottom=0.30, top=0.86)
    x = np.arange(len(values))
    colors = [GROUP_COLORS[g] for g in values["group"]]

    ax.bar(x, values["r2_full_mean"], width=0.68, color=colors, edgecolor="none", zorder=2)
    ax.errorbar(x, values["r2_full_mean"], yerr=values["r2_full_sd"], fmt="none", ecolor=de.INK,
                elinewidth=0.8, capsize=2.2, capthick=0.8, zorder=3)
    ax.scatter(x, values["r2_celltype_only"], marker="D", s=16, facecolor="white", edgecolor=de.INK,
               linewidth=0.8, zorder=4)

    ax.axhline(0, color=de.INK, lw=0.6, zorder=2.5)
    ax.set_xlim(-0.6, len(values) - 0.4)
    ax.set_ylim(-0.03, 0.32)
    ax.set_yticks(np.arange(0, 0.31, 0.05))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.2f}"))
    ax.grid(axis="y", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylabel("Cross-validated R$^2$", fontsize=FS["label"])
    ax.set_xticks(x)
    ax.set_xticklabels(values["label"], rotation=35, ha="right", rotation_mode="anchor", fontsize=FS["tick"])
    ax.tick_params(axis="x", length=0, pad=3)
    ax.tick_params(axis="y", length=2.5, width=0.6, labelsize=FS["tick"], color=de.INK)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)

    handles = [Patch(facecolor=GROUP_COLORS[g], edgecolor="none", label=g) for g in GROUPS]
    handles.append(Line2D([0], [0], marker="D", ls="none", markersize=4, markerfacecolor="white",
                          markeredgecolor=de.INK, markeredgewidth=0.8, label="Cell type only"))
    ax.legend(handles=handles, ncol=len(handles), loc="lower center", bbox_to_anchor=(0.5, 1.02), frameon=False,
              fontsize=FS["legend"], handlelength=1.1, handleheight=0.8, handletextpad=0.4, columnspacing=1.1,
              borderaxespad=0)

    de.save(fig, out_prefix)
    de.write_values(values, out_prefix)


def main() -> int:
    de.set_style()
    plot(load_values())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
