"""Paths and settings shared by the Run_PriorGRN scripts.

Motif-prior GRN methods run on the same cells and genes as the Run_BEELINE P < 0.05 runs (no expression filter,
uncorrected BEELINE GAM P < 0.05, significant TFs + the N most-varying significant genes, Fig5 TFs forced):
datasets <context>-r9-sc-cpfilt-p05-TFs{1000,2000}, contexts K562 (r9 K562 + TNF-alpha) and THP1 (THP1 + LPS).

No ATAC-seq exists for these experiments, so every method uses a motif prior without accessibility data:
  PYSCENIC     pySCENIC 0.12.1: GRNBoost2 (TF regulators) -> cisTarget pruning with the hg38 v10 gene-based motif
               rankings (500 bp up / 100 bp down and 10 kb up / 10 kb down of the TSS)
  CELLORACLE   CellOracle 0.18.0: hg38 promoter base GRN (gimmemotifs v5 scan, fpr 2%) -> bagging ridge regression
  INFERELATOR  Inferelator 0.6.4: the same CellOracle base GRN as the prior -> TFA + BBSR regression
Expression is r9's log1p layer (= log1p of its cpm layer), the matrix the BEELINE runs use.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIGS = ROOT.parent  # Fig5_GRN (since 2026-11-16 Run_PriorGRN and Run_BEELINE live inside it)
INPUT_DIR = ROOT / "input"
OUTPUT_DIR = ROOT / "output"
RES = ROOT / "resources"
INPUTS_DIR = OUTPUT_DIR / "inputs"  # per-dataset method inputs (step 01)
NETWORK_DIR = OUTPUT_DIR / "networks"  # <dataset>/<METHOD>/rankedEdges.csv, BEELINE layout and format

ADATAS = {  # symlinks to the r9 prefiltered AnnData (same files as Run_BEELINE/input)
    "K562": INPUT_DIR / "adata_prefiltered.1e4-variant.fixed.h5ad",
    "THP1": INPUT_DIR / "adata_thp1_prefiltered.fixed.h5ad",
}
BEELINE_INPUTS = FIGS / "Run_BEELINE" / "beeline" / "BEELINE-data" / "inputs" / "scRNA-Seq"
DATASETS = {ctx: [f"{ctx}-r9-sc-cpfilt-p05-TFs{n}" for n in (1000, 2000)] for ctx in ADATAS}
METHODS = ["PYSCENIC", "CELLORACLE", "INFERELATOR"]

BASE_GRN = INPUT_DIR / "hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet"  # CellOracle repo copy
CISTARGET = RES / "cistarget"
CISTARGET_DBS = [CISTARGET / f"hg38_{w}_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
                 for w in ("500bp_up_100bp_down", "10kbp_up_10kbp_down")]
MOTIF_ANNOT = CISTARGET / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
PYSCENIC_TFS = CISTARGET / "hs_hgnc_tfs.txt"

# Every TF scored by Fig5_GRN (as Run_BEELINE FORCED_TFS)
FIG5_TFS = ["RELA", "NFKB1", "JUN", "ATF2", "ELK1", "JUNB", "EGR1", "CREB1", "ATF1", "FOS", "REL", "IRF3"]
SEED = 777


def context_of(dataset: str) -> str:
    return dataset.split("-", 1)[0]


def write_ranked_edges(edges, path: Path) -> None:
    """BEELINE rankedEdges.csv: tab-separated Gene1 / Gene2 / EdgeWeight, sorted by weight (descending)."""
    e = edges[["Gene1", "Gene2", "EdgeWeight"]]
    assert (e.Gene1 != e.Gene2).all(), "self-edges"
    assert not e.duplicated(["Gene1", "Gene2"]).any(), "duplicate edges"
    e = e.sort_values("EdgeWeight", ascending=False, kind="mergesort")
    path.parent.mkdir(parents=True, exist_ok=True)
    e.to_csv(path, sep="\t", index=False)
