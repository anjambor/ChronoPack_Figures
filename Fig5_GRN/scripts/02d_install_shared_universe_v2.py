#!/usr/bin/env python
"""Check the FIG5_UNIVERSE=shared build and move it next to the with_prior outputs as the "v2" figure set.

Run after steps 02, 02b and 02c with FIG5_METHOD_SET=with_prior FIG5_UNIVERSE=shared, which write a full benchmark
tree to output/with_prior/shared_universe_build/benchmark_outputs/. This step

- checks that no Perturb-seq method (Replogle, Yao) is left in any table, that every per-TF row was scored on the
  shared universe (N = number of shared genes) and that the reference positives lie inside it;
- moves figures/5{b,c,d,e} to output/with_prior/benchmark_outputs/figures/5{b,c,d,e}_v2 and the tables to
  output/with_prior/benchmark_outputs/tables/shared_universe_v2/ (replacing earlier v2 outputs), writes the gene list
  (shared_universe_genes.txt) and a summary, and removes the build tree.

Rerunning step 02 for with_prior (full universe) clears output/with_prior/benchmark_outputs/, so run_fig5_v4.sh
rebuilds the shared set after it.
"""
from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
for p in (SCRIPTS, SCRIPTS / "capsule_lib" / "01_generate_predictions_and_metrics"):
    sys.path.insert(0, str(p))

if os.environ.get("FIG5_METHOD_SET") != "with_prior" or os.environ.get("FIG5_UNIVERSE") != "shared":
    sys.exit("run with FIG5_METHOD_SET=with_prior FIG5_UNIVERSE=shared")

import fig5_paths  # noqa: E402

import pandas as pd  # noqa: E402

import benchmark_grn_systems as b  # noqa: E402

BUILD = fig5_paths.BENCH_ROOT / "benchmark_outputs"
DEST = fig5_paths.SHARED_INSTALL_ROOT
PANELS = ["5b", "5c", "5d", "5e"]
TABLE_DEST = DEST / "tables" / "shared_universe_v2"


def check(universe: set[str]) -> dict:
    tables = sorted((BUILD / "tables").rglob("*.csv"))
    for path in tables:
        cols = pd.read_csv(path, nrows=0).columns
        if "method" in cols:
            methods = set(pd.read_csv(path, usecols=["method"])["method"].dropna().astype(str))
            bad = methods & set(fig5_paths.PERTURB_METHODS)
            assert not bad, f"{path.name}: Perturb-seq methods left: {sorted(bad)}"
        if "target" in cols:
            targets = set(pd.read_csv(path, usecols=["target"])["target"].dropna().astype(str))
            assert targets <= universe, f"{path.name}: {len(targets - universe)} targets outside the shared universe"

    per_tf = pd.read_csv(BUILD / "tables" / "publication_figures" / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv")
    scored = per_tf[per_tf["k_pos"] > 0]
    assert (scored["N"] == len(universe)).all(), "per-TF rows not scored on the shared universe"
    gts = b.load_ground_truths(universe)
    for name, gt in gts.items():
        assert set(gt["target"]) <= universe, f"{name} positives outside the universe"
    return {
        "n_tables_checked": len(tables),
        "methods": [m for m in b.METHOD_ORDER if m in set(per_tf["method"])],
        "per_tf_rows": int(len(per_tf)),
        "per_tf_rows_scored": int(len(scored)),
    }


def main() -> None:
    universe = b.load_valid_targets()  # restricted by fig5_paths.restrict_universe
    summary = {
        "universe": "shared",
        "n_genes": len(universe),
        "definition": ("BEELINE / motif-prior expression input (Run_BEELINE K562-r9-sc-cpfilt-p05-TFs2000, 2,331 genes) "
                       "that ChronoPack_K562 and ChronoPack_THP1 also rank, within Fig5's protein-coding filter"),
        "removed_methods": fig5_paths.PERTURB_METHODS,
        **check(universe),
    }
    for panel in PANELS:
        dest = DEST / "figures" / f"{panel}_v2"
        if dest.exists():
            shutil.rmtree(dest)
        shutil.move(str(BUILD / "figures" / panel), dest)
    left = sorted(p.name for p in (BUILD / "figures").iterdir())
    assert not left, f"unexpected figure folders in the build: {left}"
    if TABLE_DEST.exists():
        shutil.rmtree(TABLE_DEST)
    shutil.move(str(BUILD / "tables"), TABLE_DEST)
    (TABLE_DEST / "shared_universe_genes.txt").write_text("\n".join(sorted(universe)) + "\n")
    (TABLE_DEST / "shared_universe_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    shutil.rmtree(fig5_paths.BENCH_ROOT)
    print(json.dumps(summary, indent=2))
    print(f"Installed figures {', '.join(p + '_v2' for p in PANELS)} in {DEST / 'figures'} and tables in {TABLE_DEST}")
    print("SHARED_UNIVERSE_V2_OK")


if __name__ == "__main__":
    main()
