#!/usr/bin/env python
"""Shared paths, comparator run sets and the tie-aware ranking helper for FigS_C2C_ligand_specificity."""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "input"
FIG4_DIR = INPUT_DIR / "fig4"  # copied from ../Fig4_C2C/output (2026-10-21)
CAPSULE_DIR = INPUT_DIR / "capsule"  # the Sep 28 capsule figure's table, image and script
OUTPUT_DIR = ROOT / "output"
TABLE_DIR = OUTPUT_DIR / "tables"
FIGURE_DIR = OUTPUT_DIR / "figures"

os.environ.setdefault("MPLCONFIGDIR", str(OUTPUT_DIR / ".mplconfig"))

CHRONOPACK_LIGANDS = FIG4_DIR / "chronopack_c2c.ligand_scores.csv"
COMPARATOR_TABLES = [FIG4_DIR / f"harmonized_{n}.csv"
                     for n in ["liana_results", "liana_individual_methods", "nichenet_results", "trasig_results"]]
CAPSULE_TABLE = CAPSULE_DIR / "chronopack_specific_ligand_recovery.csv"

MISSING_PCT = 1.05  # percentile given to a ligand that no comparator run scores (as in the capsule figure)
N_LABEL = 14  # ligands labelled / shown in the bar chart (as in the capsule figure)

# Comparator runs: the 4 runs drawn in Fig. 4 (Fig4_C2C/scripts/fig4_common.method_rows).
FIG4_RUNS = {("LIANA", "static_post_baseline"), ("LIANA_CellChat", "static_post_baseline"),
             ("NicheNet", "lagged_mid_late"), ("TraSig", "official_experimental_time_pseudotime")}


def tie_ranks(df: pd.DataFrame, key: str, name_col: str) -> pd.DataFrame:
    """Rank rows by `key` (descending) with explicit ties (as C2C_information_baseline/scripts/common.tie_ranks).

    Adds rank_min (competition rank), rank_mean (expected rank under random tie breaking), rank_max and tie_size.
    Scores are rounded to 12 significant digits so floating-point noise does not split ties.
    """
    out = df.copy()
    out["_k"] = out[key].astype(float).map(lambda v: float(f"{v:.12g}") if np.isfinite(v) else -np.inf)
    out = out.sort_values(["_k", name_col], ascending=[False, True], kind="mergesort").reset_index(drop=True)
    out["position"] = np.arange(1, len(out) + 1)
    grp = out.groupby("_k", sort=False)["position"]
    out["rank_min"] = grp.transform("min")
    out["rank_max"] = grp.transform("max")
    out["rank_mean"] = grp.transform("mean")
    out["tie_size"] = out["rank_max"] - out["rank_min"] + 1
    return out.drop(columns="_k")
