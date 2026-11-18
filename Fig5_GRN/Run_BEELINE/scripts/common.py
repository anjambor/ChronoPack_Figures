"""Paths and settings shared by the Run_BEELINE scripts.

BEELINE_CONTEXT (environment, default K562) selects the cell population: K562 (r9 K562 + TNF-alpha) or THP1 (the THP1
cells of the r9 THP1 + K562 + LPS experiment).

BEELINE_GENE_FILTER (environment, default beeline) selects the expression filter applied to the r9 prefiltered genes:
  beeline     genes with counts > 0 in >= 10% of cells (BEELINE paper). Prep files in output/prep/<context>/,
              datasets <context>-r9-sc-TFs<N>.
  chronopack  no further filter: every r9 prefiltered gene, i.e. ChronoPack's own expression filter (the genes of its
              tradeSeq fit and GRN inference: detected in >= 3 cells, >= 10 total counts, technical and time-tag-bias
              genes removed; 03_PrepareForTradeSeq). Prep files in output/prep/<context>-cpfilt/, datasets
              <context>-r9-sc-cpfilt-TFs<N>.

BEELINE_PVAL (environment, default 0) is generateExpInputs.py's -p (no Bonferroni -c): 0 = no time-test filter (all TFs +
the N most-varying genes); > 0 = the paper's selection with an uncorrected GAM P cut-off (significant TFs + the N
most-varying significant genes). A nonzero cut-off adds -p<digits> to the dataset name, e.g. P = 0.05 ->
<context>-r9-sc-cpfilt-p05-TFs<N> (same prep files as P = 0).
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "input"
OUTPUT_DIR = ROOT / "output"

ADATAS = {  # symlinks to the r9 prefiltered AnnData of each context
    "K562": INPUT_DIR / "adata_prefiltered.1e4-variant.fixed.h5ad",  # 6,857 K562 + TNF-alpha cells
    "THP1": INPUT_DIR / "adata_thp1_prefiltered.fixed.h5ad",  # 787 THP1 cells of the THP1 + K562 + LPS experiment
}
CONTEXT = os.environ.get("BEELINE_CONTEXT", "K562")
if CONTEXT not in ADATAS:
    raise ValueError(f"BEELINE_CONTEXT must be one of {sorted(ADATAS)}, not {CONTEXT!r}")
ADATA = ADATAS[CONTEXT]
GENE_FILTERS = {"beeline": 0.10, "chronopack": 0.0}  # minimum fraction of cells with counts > 0
GENE_FILTER = os.environ.get("BEELINE_GENE_FILTER", "beeline")
if GENE_FILTER not in GENE_FILTERS:
    raise ValueError(f"BEELINE_GENE_FILTER must be one of {sorted(GENE_FILTERS)}, not {GENE_FILTER!r}")
FILTER_TAG = "" if GENE_FILTER == "beeline" else "-cpfilt"
PREP_DIR = OUTPUT_DIR / "prep" / f"{CONTEXT}{FILTER_TAG}"
HUMAN_TFS = INPUT_DIR / "human-tfs.csv"  # BEELINE-Networks TF list (BEELINE paper)
CHIP_GT = INPUT_DIR / "tf_gene.1kbp_200bp.stringent.tsv"  # K562 ENCODE ChIP (Fig5 ground truth), for refNetwork.csv

# BEELINE root (copied code in beeline/Beeline; data under beeline/BEELINE-data, beeline/inputs, beeline/outputs)
BEELINE_ROOT = ROOT / "beeline"
BEELINE_CODE = BEELINE_ROOT / "Beeline"
DATASET_DIR = "scRNA-Seq"

MIN_CELL_FRACTION = GENE_FILTERS[GENE_FILTER]  # beeline: drop genes expressed in fewer than 10% of cells
# BEELINE's Bonferroni P < 0.01 time test keeps only 277 of 3,611 genes on the r9 K562 single cells, too few for TFs + N,
# so the p-value filter is off (generateExpInputs.py -p 0): all TFs + the N most-varying non-TF genes (2026-10-19).
# 2026-10-21: P = 0.05 (uncorrected) runs added for the no-filter (chronopack) gene set.
PVAL = float(os.environ.get("BEELINE_PVAL", "0"))
PVAL_TAG = "" if PVAL == 0 else "-p" + f"{PVAL:g}".split(".")[-1]
N_GENES = [1000, 2000]  # TFs + N
# Every TF scored by Fig5_GRN (CORE_TFS, K562/THP1 evaluation TFs); added even if not selected (both contexts)
FORCED_TFS = ["RELA", "NFKB1", "JUN", "ATF2", "ELK1", "JUNB", "EGR1", "CREB1", "ATF1", "FOS", "REL", "IRF3"]


def dataset_name(n: int, context: str = CONTEXT, filter_tag: str = FILTER_TAG, pval_tag: str = PVAL_TAG) -> str:
    return f"{context}-r9-sc{filter_tag}{pval_tag}-TFs{n}"
