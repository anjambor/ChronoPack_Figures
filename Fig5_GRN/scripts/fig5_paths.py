"""Fig5_GRN (v4; folder renamed from Fig5_GRN_v4 on 2026-11-17) directory layout, shared by the Fig5 scripts and the copied capsule scripts.

The scripts under capsule_lib/ are copies of grn_inference_v2/reproduction/scripts with
their path roots redirected here: inputs come from Fig5_GRN/input/ and all generated
tables/figures go to Fig5_GRN_v4/output/benchmark_outputs/ (same sub-layout as the capsule).
v3 differs from ../Fig5_GRN_v2 in the BEELINE predictions (BEELINE_RESULTS_DIR: the P < 0.05 run) and adds
the motif-prior methods (PRIOR_RESULTS_DIR) in the "with_prior" method set.
v4 differs from ../Fig5_GRN_v3 only in ChronoPack's TF rankings: step 01 computes them with the current ChronoPack
(cp.tl.grn.rank_regulators defaults since commit e0c89fb, 2026-10-22) from the pinned OmniPath snapshot
(OMNIPATH_SNAPSHOT) instead of reading the 2026-09-30 r9 rankings. A requested TF that a method has no predictions
for is left out of that method's scores, as in v3; with the new rankings this now applies to ChronoPack_THP1 CREB1
(priority_act 0 for TLR4, so infer drops it). CREB1 is still scored for every other method.

FIG5_METHOD_SET selects the comparator methods of every cross-method panel:
  expression  (default) ChronoPack, BEELINE (GENIE3, SINCERITIES, GRNBOOST2, PPCOR), Replogle, Yao: v2's methods.
              Benchmark outputs in output/benchmark_outputs/.
  with_prior  the same plus the K562 motif-prior methods (pySCENIC, CellOracle, Inferelator; Run_PriorGRN).
              Benchmark outputs in output/with_prior/benchmark_outputs/.
ChronoPack's own predictions (step 01, output/chronopack/) are shared by both sets.

FIG5_UNIVERSE selects the evaluation gene universe (user request 2026-11-11; with_prior only):
  full    (default) all filtered protein-coding genes (19,338 targets), every method incl. Replogle and Yao.
  shared  only the genes every non-Perturb-seq method can rank: the BEELINE / motif-prior expression input
          (2,331 genes) that ChronoPack_K562 and ChronoPack_THP1 also rank, within the protein-coding filter;
          Replogle and Yao are removed from every panel. Built in output/with_prior/shared_universe_build/ and
          moved by step 02d to output/with_prior/benchmark_outputs/figures/5{b,c,d,e}_v2 and tables/shared_universe_v2.
"""
from __future__ import annotations

import os
from pathlib import Path

FIG_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = FIG_DIR / "input"
CAPSULE_INPUT_DIR = INPUT_DIR / "capsule"  # files copied from the GRN reproduction capsule
EXTERNAL_DIR = CAPSULE_INPUT_DIR / "external"  # mirrors /new-stg/home/alex/{ref,downloads,projects}
# BEELINE GENIE3 / GRNBOOST2 / PPCOR / SINCERITIES predictions: the r9 single-cell P < 0.05 TFs + 2000 run of
# Run_BEELINE (symlink; no expression filter, uncorrected GAM P < 0.05, significant TFs + 2,000 most-varying
# significant genes, Fig5 TFs forced), replacing v2's 10%-filter run (K562-r9-sc-TFs2000).
BEELINE_RESULTS_DIR = INPUT_DIR / "beeline_r9" / "K562-r9-sc-cpfilt-p05-TFs2000"
# pySCENIC / CellOracle / Inferelator networks on the same dataset (Run_PriorGRN, symlink)
PRIOR_RESULTS_DIR = INPUT_DIR / "prior_grn" / "K562-r9-sc-cpfilt-p05-TFs2000"
PRIOR_METHODS = ["PYSCENIC", "CELLORACLE", "INFERELATOR"]
METHOD_SET = os.environ.get("FIG5_METHOD_SET", "expression")
if METHOD_SET not in ("expression", "with_prior"):
    raise ValueError(f"FIG5_METHOD_SET must be expression or with_prior, not {METHOD_SET!r}")
WITH_PRIOR = METHOD_SET == "with_prior"
CHRONOPACK_INPUT_DIR = INPUT_DIR / "chronopack"  # symlinks to ChronoPack r9 data / resources
# OmniPath snapshot the ChronoPack analyses use (ChronoPack/analyses/resources/, symlink): rankings need no network
OMNIPATH_SNAPSHOT = CHRONOPACK_INPUT_DIR / "omnipath_snapshot_2025-10-21"

# FIG5_OUTPUT_DIR redirects all outputs (default: Fig5_GRN/output).
OUTPUT_DIR = Path(os.environ.get("FIG5_OUTPUT_DIR", FIG_DIR / "output")).resolve()
CHRONOPACK_OUTPUT_DIR = OUTPUT_DIR / "chronopack"
UNIVERSE = os.environ.get("FIG5_UNIVERSE", "full")
if UNIVERSE not in ("full", "shared"):
    raise ValueError(f"FIG5_UNIVERSE must be full or shared, not {UNIVERSE!r}")
SHARED_UNIVERSE = UNIVERSE == "shared"
if SHARED_UNIVERSE and not WITH_PRIOR:
    raise ValueError("FIG5_UNIVERSE=shared is only set up for FIG5_METHOD_SET=with_prior")
PERTURB_METHODS = ["Replogle2022_K562", "Yao2023_THP1"]
BENCH_ROOT = OUTPUT_DIR / "with_prior" if WITH_PRIOR else OUTPUT_DIR  # capsule scripts write to BENCH_ROOT / "benchmark_outputs" / ...
SHARED_INSTALL_ROOT = BENCH_ROOT / "benchmark_outputs"  # step 02d moves the shared-universe outputs here
if SHARED_UNIVERSE:
    BENCH_ROOT = BENCH_ROOT / "shared_universe_build"
TABLE_DIR = BENCH_ROOT / "benchmark_outputs" / "tables"
PUBLICATION_TABLE_DIR = TABLE_DIR / "publication_figures"
FIGURE_DIR = BENCH_ROOT / "benchmark_outputs" / "figures"

CONTEXTS = {
    # context: (method name, OmniPath source, is_ligand)
    "K562": ("ChronoPack_K562", "TNF", True),
    "THP1": ("ChronoPack_THP1", "TLR4", False),
}

# Motif-prior methods (FIG5_METHOD_SET=with_prior): display names and colours (purple/brown, distinct from v2's
# ChronoPack orange, BEELINE greens, Replogle pink and Yao blue)
PRIOR_LABELS = {"PYSCENIC": "pySCENIC", "CELLORACLE": "CellOracle", "INFERELATOR": "Inferelator"}
PRIOR_COLORS = {"PYSCENIC": "#6A3D9A", "CELLORACLE": "#9E7CC1", "INFERELATOR": "#8C564B"}


def insert_prior(methods, after: str = "PPCOR") -> list:
    """`methods` with PRIOR_METHODS inserted after `after` (the last BEELINE method) when FIG5_METHOD_SET=with_prior,
    and without the Perturb-seq methods when FIG5_UNIVERSE=shared."""
    methods = list(methods)
    if WITH_PRIOR and after in methods:
        i = methods.index(after) + 1
        methods[i:i] = [m for m in PRIOR_METHODS if m not in methods]
    if SHARED_UNIVERSE:
        methods = [m for m in methods if m not in PERTURB_METHODS]
    return methods


def shared_universe_genes() -> set[str]:
    """Genes every non-Perturb-seq method can rank: the BEELINE expression input (identical to the motif-prior
    methods' input) that ChronoPack_K562 and ChronoPack_THP1 also rank (step 01 targets)."""
    global _SHARED_GENES
    if _SHARED_GENES is None:
        import pandas as pd

        ds = BEELINE_RESULTS_DIR.name
        beeline_expr = BEELINE_RESULTS_DIR.resolve().parents[2] / "inputs" / "scRNA-Seq" / ds / "GENIE3" / "ExpressionData.csv"
        prior_expr = PRIOR_RESULTS_DIR.resolve().parents[1] / "inputs" / ds / "expr_log1p.tsv"
        genes = set(pd.read_csv(beeline_expr, sep="\t", nrows=0, index_col=0).columns)
        if set(pd.read_csv(prior_expr, sep="\t", nrows=0, index_col=0).columns) != genes:
            raise ValueError("BEELINE and motif-prior expression inputs have different genes")
        for context in CONTEXTS:
            genes &= set(pd.read_csv(CHRONOPACK_OUTPUT_DIR / context / "grn.edges.csv", usecols=["target"])["target"].astype(str))
        _SHARED_GENES = genes
    return set(_SHARED_GENES)


_SHARED_GENES = None


def restrict_universe(valid_targets: set[str]) -> set[str]:
    """The evaluation universe: `valid_targets` (protein-coding filter), cut to shared_universe_genes() when
    FIG5_UNIVERSE=shared."""
    return valid_targets & shared_universe_genes() if SHARED_UNIVERSE else valid_targets
