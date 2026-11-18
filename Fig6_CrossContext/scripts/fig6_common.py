#!/usr/bin/env python3
"""Shared constants and loaders for the Fig6 (originally Fig7) ORA + covariate + PLS motif dot plot.

Adapted from
/new-stg/home/alex/tmp/fig3_dm_dotplot/ora_v3/analysis/plot_combined_ora_pls_motif_dotplot.py.
The motif panel is rebuilt directly from the raw motif summary and TATA positional
summary, porting the logic of
/new-stg/home/alex/tmp/motif_enrichment/input_S8/create_s8_crematched_tata_pls_specific_strict_horizontal_pls_only.py
(and the modules it loads) instead of reading its pre-built .data.csv/.rows.csv.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "input"
OUTPUT_DIR = ROOT / "output"
ORA_INPUT = INPUT_DIR / "ORA_results_v2_newbg_newfg4.csv"
ORA_CURATED = INPUT_DIR / "curated_ordered_grouped_specific_terms.csv"
MOTIF_INPUT = INPUT_DIR / "motif_summary_S10_CRE-matched_with_counts.csv"
MOTIF_VARIANT = "S10_CRE-matched"
TATA_INPUT = INPUT_DIR / "tata_box_positional_S10_summary_with_counterfactual_q.csv"
TATA_BACKGROUND = "PLS-specific"
TATA_FILTER = "strict"
TATA_Q_COLUMN = "median_bh_q_if_positional"

MODE_LEVELS = [
    "impulse-up_1",
    "impulse-up_2",
    "trend-up",
    "trend-down_1",
    "trend-down_2",
    "transient-down",
]
CELL_LEVELS = ["K562a", "THP1"]
MODE_LABELS = {
    "impulse-up_1": "Impulse Up 1",
    "impulse-up_2": "Impulse Up 2",
    "trend-up": "Trend Up",
    "trend-down_1": "Trend Down 1",
    "trend-down_2": "Trend Down 2",
    "transient-down": "Transient Down",
}
ORA_Q_CUTOFF = 0.05
MOTIF_Q_CUTOFF = 0.10
COVARIATE_PADJ_CUTOFF = 0.10
EXCLUDE_ORA_DESCRIPTIONS = {
    "Tullai_IEG",
    "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
    "GOBP_CHROMATIN_REMODELING",
}

# Row definitions from the reference ...covariates_PLS_only.rows.csv, in plot order.
COVARIATE_ROWS = [
    ("gene_length_kb", "Gene length"),
    ("median_log10_decay_rate", "mRNA decay rate"),
    ("is_ARE", "ARE motif"),
]
PLS_MOTIF_ROWS = [
    ("jaspar", "MA0107.1_RELA", "NF-kB (RELA)"),
    ("homerdefault", "AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer", "TRE/AP-1 (AP-1)"),
    ("homerdefault", "CRE(bZIP)/Promoter/Homer", "CRE"),
    ("homerdefault", "TATA-Box(TBP)/Promoter/Homer", "TATA box"),
]
TATA_DISPLAY = "TATA box"


def neg_log10(pvalue: float | None) -> float:
    if pvalue is None or pd.isna(pvalue) or pvalue <= 0:
        return np.nan
    return -math.log10(float(pvalue))


def motif_area(nlp: float | None) -> float:
    if nlp is None or pd.isna(nlp):
        return 0.0
    capped = min(float(nlp), 3.0)
    return 24.0 + 95.0 * capped


def make_row_geometry() -> tuple[dict[tuple[str, str], float], list[tuple[str, str]], dict[str, float]]:
    contexts = [(mode, cell) for mode in MODE_LEVELS for cell in CELL_LEVELS]
    y_lookup = {}
    for i, (mode, cell) in enumerate(contexts):
        y_lookup[(cell, mode)] = len(contexts) - 1 - i
    mode_mid = {
        mode: (y_lookup[("K562a", mode)] + y_lookup[("THP1", mode)]) / 2
        for mode in MODE_LEVELS
    }
    return y_lookup, contexts, mode_mid


def load_ora_terms() -> pd.DataFrame:
    curated = pd.read_csv(ORA_CURATED)
    curated = curated[~curated["Description"].isin(EXCLUDE_ORA_DESCRIPTIONS)].copy()
    curated["panel"] = "ORA"
    curated["input_order"] = np.arange(len(curated))
    return curated


def load_ora_plot(ora_terms: pd.DataFrame, y_lookup: dict[tuple[str, str], float]) -> pd.DataFrame:
    df = pd.read_csv(ORA_INPUT)
    rows = []
    for _, term in ora_terms.iterrows():
        for mode in MODE_LEVELS:
            for cell in CELL_LEVELS:
                hit = df[
                    (df["Description"] == term["Description"])
                    & (df["cluster"] == mode)
                    & (df["cell_type"] == cell)
                ]
                rec = {
                    "display": term["display"],
                    "group": term["mode"],
                    "cell": cell,
                    "cluster": mode,
                    "y": y_lookup[(cell, mode)],
                    "present": False,
                    "GeneRatio_val": np.nan,
                    "FoldEnrichment": np.nan,
                    "log2FE": np.nan,
                    "qvalue": np.nan,
                }
                if not hit.empty:
                    row = hit.iloc[0]
                    fe = float(row["FoldEnrichment"])
                    rec.update(
                        present=True,
                        GeneRatio_val=float(row["GeneRatio_val"]),
                        FoldEnrichment=fe,
                        log2FE=math.log2(fe) if fe > 0 else np.nan,
                        qvalue=float(row["qvalue"]),
                    )
                rows.append(rec)
    return pd.DataFrame(rows)


def load_motif_rows() -> pd.DataFrame:
    records = [
        {"group": "Covariates", "cre_class": "Covariate", "predictor": predictor, "display": label}
        for predictor, label in COVARIATE_ROWS
    ]
    records += [
        {"group": "PLS motifs", "cre_class": "PLS", "motif_db": db, "motif_name": name, "display": label}
        for db, name, label in PLS_MOTIF_ROWS
    ]
    return pd.DataFrame(records)


def load_pls_motif_plot(
    y_lookup: dict[tuple[str, str], float],
    motif_csv: Path = MOTIF_INPUT,
    variant: str = MOTIF_VARIANT,
    tata_csv: Path = TATA_INPUT,
    tata_q_column: str | None = TATA_Q_COLUMN,
    q_cutoff: float = MOTIF_Q_CUTOFF,
) -> pd.DataFrame:
    """PLS motif points: log2 Fisher enrichment, -log10 Fisher p, outline on two-sided q.

    The TATA box row is replaced by the positional TATA summary (PLS-specific, strict).
    Its outline uses ``tata_q_column`` from that file; with ``None`` it keeps the HOMER
    TATA-Box q, as the reference S8 builder did.
    """
    specs = pd.DataFrame(PLS_MOTIF_ROWS, columns=["motif_db", "motif_name", "display"])
    motif = pd.read_csv(motif_csv)
    motif = motif[
        motif["variant"].eq(variant)
        & motif["cre_class"].eq("PLS")
        & motif["cell"].isin(CELL_LEVELS)
        & motif["cluster"].isin(MODE_LEVELS)
    ]
    motif = motif.merge(specs, on=["motif_db", "motif_name"], how="inner")
    motif = motif[
        ["cell", "cluster", "motif_db", "motif_name", "display",
         "median_log2_enrichment_fisher", "median_fisher_p", "median_fisher_two_sided_q"]
    ].copy()
    motif["effect"] = motif["median_log2_enrichment_fisher"]
    motif["p_value"] = motif["median_fisher_p"]
    motif["q_value"] = motif["median_fisher_two_sided_q"]

    tata = pd.read_csv(tata_csv)
    tata = tata[tata["background"].eq(TATA_BACKGROUND) & tata["filter"].eq(TATA_FILTER)]
    tata_cols = ["cell", "cluster", "median_log2_enr", "median_fisher_p"]
    if tata_q_column is not None:
        tata_cols.append(tata_q_column)
    motif = motif.merge(
        tata[tata_cols].rename(columns={"median_fisher_p": "tata_median_fisher_p"}),
        on=["cell", "cluster"],
        how="left",
    )
    is_tata = motif["display"].eq(TATA_DISPLAY)
    motif.loc[is_tata, "effect"] = motif.loc[is_tata, "median_log2_enr"]
    motif.loc[is_tata, "p_value"] = motif.loc[is_tata, "tata_median_fisher_p"]
    if tata_q_column is not None:
        motif.loc[is_tata, "q_value"] = motif.loc[is_tata, tata_q_column]
    motif["tata_variant_applied"] = is_tata

    motif["plot_neg_log10_p"] = motif["p_value"].map(neg_log10)
    motif["outlined"] = motif["q_value"].notna() & (motif["q_value"] < q_cutoff)
    motif["area"] = motif["plot_neg_log10_p"].map(motif_area)
    motif["y"] = motif.apply(lambda r: y_lookup[(r["cell"], r["cluster"])], axis=1)
    if len(motif) != len(PLS_MOTIF_ROWS) * len(CELL_LEVELS) * len(MODE_LEVELS):
        raise ValueError(f"Expected one row per motif/cell/mode, got {len(motif)}")
    return motif.drop(columns=["median_log2_enr", "tata_median_fisher_p"])


def load_covariate_plot(
    y_lookup: dict[tuple[str, str], float],
    covariate_csv: Path,
    padj_cutoff: float = COVARIATE_PADJ_CUTOFF,
) -> pd.DataFrame:
    labels = pd.DataFrame(COVARIATE_ROWS, columns=["predictor", "display"])
    cov = pd.read_csv(covariate_csv)
    cov = cov[
        cov["context"].isin(CELL_LEVELS)
        & cov["foreground"].isin(MODE_LEVELS)
        & cov["predictor"].isin(labels["predictor"])
    ].copy()
    cov = cov.merge(labels, on="predictor", how="inner")
    cov["cell"] = cov["context"]
    cov["cluster"] = cov["foreground"]
    cov["effect"] = cov["coef"] / math.log(2)
    cov["plot_neg_log10_p"] = cov["p_value"].map(neg_log10)
    cov["q_value"] = cov["p_adj_global"]
    cov["outlined"] = cov["p_adj_global"].notna() & (cov["p_adj_global"] < padj_cutoff)
    cov["area"] = cov["plot_neg_log10_p"].map(motif_area)
    cov["y"] = cov.apply(lambda r: y_lookup[(r["cell"], r["cluster"])], axis=1)
    if len(cov) != len(COVARIATE_ROWS) * len(CELL_LEVELS) * len(MODE_LEVELS):
        raise ValueError(f"Expected one row per covariate/cell/mode, got {len(cov)}")
    return cov


def load_motif_plot(
    y_lookup: dict[tuple[str, str], float],
    covariate_csv: Path,
    motif_q_cutoff: float = MOTIF_Q_CUTOFF,
    covariate_padj_cutoff: float = COVARIATE_PADJ_CUTOFF,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    return (
        load_pls_motif_plot(y_lookup, q_cutoff=motif_q_cutoff),
        load_covariate_plot(y_lookup, covariate_csv, padj_cutoff=covariate_padj_cutoff),
    )
