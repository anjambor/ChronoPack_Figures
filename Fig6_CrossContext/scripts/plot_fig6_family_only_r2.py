#!/usr/bin/env python3
"""Per-family cross-validated R^2: how well each feature family predicts each kinetic outcome on its own.

Each family-only elastic-net model (cv_summary_results.csv, level ``family_only``) uses one feature family plus
cell type: the v4 screen partials cell type out of the response and features, fits the elastic net on the
residuals and adds the cell-type component back for prediction (Fig6_Reproduction/docs/methods.md). Its R^2
therefore includes what cell type alone explains, which dominates the timing outcomes (t_peak: cell type alone
0.28; every family-only model 0.28-0.30).

- ``fig6_family_only_cv_r2_gain`` (main figure): family-only mean CV R^2 minus the cell-type-only mean CV R^2
  (cv_celltype_only_results.csv, held-out intercept + cell-type fits on the same outcome-specific folds), i.e. the
  R^2 a family adds over cell type when used alone. Small negative values mean the family's features cost more
  in fitting noise than they add.
- ``fig6_family_only_cv_r2_raw``: the family-only mean CV R^2 itself, with the cell-type-only model as the
  top row for reference.

Same families, row order, outcome columns and heatmap style as 6d (plot_fig6de.py), so the two read side by
side: 6d = signal a family carries that the others do not; this figure = signal a family carries on its own.
Writes output/fig6_family_only_cv_r2_{gain,raw}.{pdf,png} and a shared _values.csv.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "output" / ".mplconfig"))

import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fig6_common as combined  # noqa: E402
import plot_fig6de as de  # noqa: E402
import plot_fig6_full_model_r2 as full  # noqa: E402

OUT_GAIN = combined.OUTPUT_DIR / "fig6_family_only_cv_r2_gain"
OUT_RAW = combined.OUTPUT_DIR / "fig6_family_only_cv_r2_raw"
VALUES = combined.OUTPUT_DIR / "fig6_family_only_cv_r2_values.csv"
CELLTYPE_ROW = "Cell type only"

GAIN_LABEL = "CV R² gain over cell type alone (family-only model)"
GAIN_LABEL_WRAPPED = "CV R² gain over\ncell type alone"
RAW_LABEL = "CV R² of the family-only model"
RAW_LABEL_WRAPPED = "CV R² of the\nfamily-only model"


def load_values() -> pd.DataFrame:
    cv = pd.read_csv(de.CV_INPUT)
    de.check_family_labels(cv)
    fam = cv[cv["level"].eq("family_only")].set_index(["family", "outcome"])
    ct = pd.read_csv(full.CELLTYPE_INPUT).set_index("outcome")
    rows = []
    for family, label in de.FAMILIES:
        for outcome in de.OUTCOMES:
            if (family, outcome) not in fam.index or outcome not in ct.index:
                raise ValueError(f"missing family_only or cell-type-only value for {family} / {outcome}")
            r = fam.loc[(family, outcome)]
            rows.append({
                "family": family, "family_label": label, "outcome": outcome,
                "r2_family_only_mean": r["r2_mean"], "r2_family_only_sd": r["r2_sd"],
                "r2_family_only_min": r["r2_min"], "r2_family_only_max": r["r2_max"], "n_repeats": r["n_repeats"],
                "r2_celltype_only": ct.loc[outcome, "r2_mean"], "r2_celltype_only_sd": ct.loc[outcome, "r2_sd"],
            })
    out = pd.DataFrame(rows)
    out["r2_gain_over_celltype"] = out["r2_family_only_mean"] - out["r2_celltype_only"]
    return out


def plot(mat: pd.DataFrame, row_labels: list[str], label: str, label_wrapped: str, out_prefix: Path) -> None:
    """Standalone heatmap with 6d's geometry (rows scale the figure height so cells keep 6d's size)."""
    n_extra = len(row_labels) - len(de.FAMILIES)
    base_h, row_h = 3.5, 0.40
    fig_h = base_h + n_extra * row_h
    fig, ax = plt.subplots(figsize=(8.4, fig_h))
    bottom, top = 0.27 * base_h / fig_h, 1 - 0.04 * base_h / fig_h
    fig.subplots_adjust(left=0.245, right=0.85, top=top, bottom=bottom)
    cax = fig.add_axes([0.867, bottom + 0.08 * base_h / fig_h, 0.013, top - bottom - 0.11 * base_h / fig_h])
    de.draw_d(ax, cax, mat, de.STANDALONE, "vertical", row_labels=row_labels, label=label, label_wrapped=label_wrapped)
    de.save(fig, out_prefix)


def main() -> int:
    de.set_style()
    values = load_values()
    gain = values.pivot(index="family", columns="outcome", values="r2_gain_over_celltype")
    gain = gain.reindex(index=[f for f, _ in de.FAMILIES], columns=de.OUTCOMES)
    plot(gain, [lab for _, lab in de.FAMILIES], GAIN_LABEL, GAIN_LABEL_WRAPPED, OUT_GAIN)

    raw = values.pivot(index="family", columns="outcome", values="r2_family_only_mean")
    raw = raw.reindex(index=[f for f, _ in de.FAMILIES], columns=de.OUTCOMES)
    ct = values.drop_duplicates("outcome").set_index("outcome")["r2_celltype_only"].reindex(de.OUTCOMES)
    raw = pd.concat([ct.to_frame(CELLTYPE_ROW).T, raw])
    plot(raw, [CELLTYPE_ROW] + [lab for _, lab in de.FAMILIES], RAW_LABEL, RAW_LABEL_WRAPPED, OUT_RAW)

    values.to_csv(VALUES, index=False)
    print(f"Wrote {VALUES}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
