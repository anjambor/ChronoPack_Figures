#!/usr/bin/env python
"""Shared paths, constants and top-10 selection for Fig4 v2.

The ligand rankings are v1's (../Fig4_C2C/output, copied to input/v1/); only the trajectory
panels change. Constants and the top-10 selection are those of Fig4_C2C/scripts/07_plot_fig4.py.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

FIG_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = FIG_DIR / "input"
V1_DIR = INPUT_DIR / "v1"
OUTPUT_DIR = FIG_DIR / "output"
CHRONOPACK_DIR = OUTPUT_DIR / "chronopack"
FIGURE_DIR = OUTPUT_DIR / "figures"

os.environ.setdefault("MPLCONFIGDIR", str(OUTPUT_DIR / ".mplconfig"))

ADATA = INPUT_DIR / "adata_full.chk1.h5ad"
THP1_FIT_DIR = INPUT_DIR / "thp1_tradeseq"  # -> r9 intmd/sc-THP1_to_tradeSeq.fixed
K562_FIT_DIR = INPUT_DIR / "k562_tradeseq"  # -> r9 intmd/sc-K562_to_tradeSeq.fixed
STORED_SCE = "sce.1trajectory.nknots8.rds"  # r9 final fit (genes with p < 0.05)
NKNOTS = 8
TIME_MAX = 165.0  # time tags 1-12 = 0-165 min, 15-min intervals
N_POINTS = 166  # 1-min grid

CURVES_CSV = CHRONOPACK_DIR / "fitted_curves.csv"
FIT_SUMMARY_CSV = CHRONOPACK_DIR / "ligand_fit_summary.csv"

SENDER = "THP1"
RECEIVER = "K562"
TIME_KEY = "time"
CELL_TYPE_KEY = "cell_type"
EXPR_LAYER = "log1p"
BASE_INTERVAL_MIN = 15
CAUSAL_WINDOW = (0.0, 45.0)
TRASIG_N_PERMS = 10000
POSITIVE_CONTROL_PAIRS = ["TNF->TNFRSF1B", "TNF->TNFRSF1A"]
RESPONSE_GENE = "NFKBIA"
TOP_N = 10

METHODS = [
    {"label": "ChronoPack", "color": "#2f6f4e"},
    {"label": "LIANA", "color": "#4c78a8"},
    {"label": "CellChat", "color": "#f28e2b"},
    {"label": "NicheNet", "color": "#c44e52"},
    {"label": "TraSig", "color": "#7f5aa2"},
]
CURVE_METHODS = [spec for spec in METHODS if spec["label"] != "CellChat"]


def time_to_minute(values):
    vals = pd.to_numeric(pd.Series(values).astype(str), errors="coerce")
    return (vals - 1) * BASE_INTERVAL_MIN


def to_dense(x):
    return x.toarray() if hasattr(x, "toarray") else np.asarray(x)


def top_unique_ligands(rows: pd.DataFrame, n: int) -> pd.DataFrame:
    rows = rows.sort_values("rank", kind="stable").copy()
    out = rows.drop_duplicates("ligand", keep="first").head(n).copy()
    out["ligand_rank"] = np.arange(1, len(out) + 1)
    return out


def chronopack_rows() -> pd.DataFrame:
    ligands = pd.read_csv(V1_DIR / "chronopack_c2c.ligand_scores.csv")
    pairs = pd.read_csv(V1_DIR / "chronopack_c2c.interactions.csv").set_index("pair")
    rows = ligands.copy()
    rows["pair"] = rows["best_pair"]
    rows["receptor"] = rows["best_pair"].map(pairs["receptor"])
    rows["activity_sum"] = rows["max_lr_activity"]
    rows["logp"] = rows["ligand_significance"]
    rows["rec_det"] = rows["best_pair"].map(pairs["rec_det"])
    rows["n_window_points"] = rows["best_pair"].map(pairs["n_window_points"])
    rows["rank"] = rows["zscore_composite_rank"].astype(int)
    rows["interaction_name"] = rows["best_pair"]
    rows["score"] = rows["composite_score"]
    rows["score_source"] = "chronopack_zscore_ligand_significance_plus_max_lr_activity"
    rows["method"] = "ChronoPack"
    rows["run_id"] = "chronopack_c2c_infer"
    return rows[["ligand", "receptor", "pair", "activity_sum", "rec_det", "logp", "n_window_points", "rank", "interaction_name", "score", "score_source", "method", "run_id"]]


def method_rows(label: str) -> pd.DataFrame:
    if label == "ChronoPack":
        return chronopack_rows()
    if label == "LIANA":
        rows = pd.read_csv(V1_DIR / "harmonized_liana_results.csv")
        return rows.loc[(rows["method"] == "LIANA") & (rows["run_id"] == "static_post_baseline")].copy()
    if label == "CellChat":
        rows = pd.read_csv(V1_DIR / "harmonized_liana_individual_methods.csv")
        return rows.loc[(rows["method"] == "LIANA_CellChat") & (rows["run_id"] == "static_post_baseline")].copy()
    if label == "NicheNet":
        return pd.read_csv(V1_DIR / "harmonized_nichenet_results.csv")
    if label == "TraSig":
        return pd.read_csv(V1_DIR / "harmonized_trasig_results.csv")
    raise ValueError(label)


def build_top_manifest() -> pd.DataFrame:
    parts = []
    for spec in METHODS:
        rows = top_unique_ligands(method_rows(spec["label"]), n=TOP_N)
        rows["plot_method"] = spec["label"]
        parts.append(rows)
    return pd.concat(parts, ignore_index=True)


def curve_ligands() -> list[str]:
    """Ligands drawn in the trajectory panels: the top 10 of every method except CellChat."""
    manifest = build_top_manifest()
    keep = manifest["plot_method"].isin([spec["label"] for spec in CURVE_METHODS])
    return sorted(set(manifest.loc[keep, "ligand"].astype(str)))
