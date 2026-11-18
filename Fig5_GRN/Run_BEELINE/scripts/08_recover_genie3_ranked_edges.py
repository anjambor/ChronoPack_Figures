#!/usr/bin/env python
"""Write rankedEdges.csv for a GENIE3 run whose network was computed but whose job hung at dask shutdown.

BEELINE's runArboreto.py (arboreto.sif) computes the whole network, writes outFile.txt with one
DataFrame.to_csv call, and only then exits; BLRunner's genie3Runner.parseOutput converts outFile.txt to
rankedEdges.csv. THP1-r9-sc-TFs2000 GENIE3 (job 941868, bigcompute) wrote outFile.txt at 19:45 on 2026-10-20, then hung
in distributed's cluster close (OSError: Timed out trying to connect to inproc://..., at 19:46), so parseOutput never
ran; the job was cancelled on 2026-10-21.

Checks on outFile.txt (any failure aborts without writing):
  header TF / target / importance; ends with a newline; TF and target names are the dataset's genes; no self-edges,
  no duplicate pairs; every gene is a target; importances finite, > 0 and sorted descending (arboreto's order); the
  number of pairs missing from all G x (G - 1) is small (arboreto drops zero-importance links; the completed
  THP1-r9-sc-TFs1000 run misses 63).
Conversion: genie3Runner.parseOutput's code, unchanged (pd.read_csv, then f"{TF}\\t{target}\\t{importance}" per row).
It is first applied to a completed run (--reference) and must reproduce that run's BEELINE-written rankedEdges.csv
byte for byte. Run with the BEELINE env python (same pandas as BLRunner):
  python scripts/08_recover_genie3_ranked_edges.py THP1-r9-sc-TFs2000 --reference THP1-r9-sc-TFs1000
"""
from __future__ import annotations

import argparse
import filecmp
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.BEELINE_ROOT / "outputs" / c.DATASET_DIR
INP = c.BEELINE_ROOT / "BEELINE-data" / "inputs" / c.DATASET_DIR


def parse_output(output_file: Path, ranked_edges_file: Path) -> None:
    """genie3Runner.parseOutput, verbatim apart from the paths."""
    OutDF = pd.read_csv(output_file, sep='\t', header=0)
    with open(ranked_edges_file, 'w') as outFile:
        outFile.write('Gene1\tGene2\tEdgeWeight\n')
        for idx, row in OutDF.iterrows():
            outFile.write(f"{row['TF']}\t{row['target']}\t{row['importance']}\n")


def check(dataset: str) -> list[str]:
    f = OUT / dataset / "GENIE3" / "outFile.txt"
    genes = pd.read_csv(INP / dataset / "ExpressionData.csv", index_col=0, usecols=[0]).index
    G = len(genes)
    with open(f, "rb") as h:
        h.seek(-1, 2)
        assert h.read(1) == b"\n", "outFile.txt does not end with a newline"
    df = pd.read_csv(f, sep="\t", header=0)
    assert list(df.columns) == ["TF", "target", "importance"], df.columns
    gs = set(genes)
    assert set(df.TF) <= gs and set(df.target) <= gs, "names outside the dataset's genes"
    assert (df.TF != df.target).all(), "self-edges"
    assert not df.duplicated(["TF", "target"]).any(), "duplicate pairs"
    assert set(df.target) == gs, f"{len(gs - set(df.target))} genes never a target"
    imp = df.importance.values.astype(float)
    assert np.isfinite(imp).all() and (imp > 0).all(), "non-finite or non-positive importances"
    assert (np.diff(imp) <= 0).all(), "importances not sorted descending"
    missing = G * (G - 1) - len(df)
    assert 0 <= missing <= 0.001 * G * (G - 1), f"{missing} pairs missing"
    per_target = df.groupby("target").size()
    return [f"{dataset}: {G} genes, {len(df):,} edges ({missing} of {G * (G - 1):,} pairs absent = zero importance), "
            f"{df.TF.nunique()} regulators, edges per target {per_target.min()}-{per_target.max()}, "
            f"importance {imp.min():.3g}-{imp.max():.3g}"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("--reference", required=True, help="completed GENIE3 run whose rankedEdges.csv BEELINE wrote")
    a = ap.parse_args()
    lines = check(a.reference) + check(a.dataset)

    ref_dir = OUT / a.reference / "GENIE3"
    with tempfile.TemporaryDirectory(dir=c.BEELINE_ROOT / "tmp") as tmp:
        t = Path(tmp) / "rankedEdges.csv"
        parse_output(ref_dir / "outFile.txt", t)
        same = filecmp.cmp(t, ref_dir / "rankedEdges.csv", shallow=False)
    lines.append(f"conversion reproduces {a.reference}/GENIE3/rankedEdges.csv byte for byte: {same}")
    assert same, "conversion does not reproduce BEELINE's rankedEdges.csv"

    out_dir = OUT / a.dataset / "GENIE3"
    target = out_dir / "rankedEdges.csv"
    assert not target.exists(), f"{target} exists; not overwriting"
    parse_output(out_dir / "outFile.txt", target)
    lines.append(f"wrote {target} ({sum(1 for _ in open(target)) - 1:,} edges)")
    note = ("rankedEdges.csv written by Run_BEELINE/scripts/08_recover_genie3_ranked_edges.py (2026-10-21), with "
            "genie3Runner.parseOutput's code, from the outFile.txt of job 941868. The job hung after writing "
            "outFile.txt (dask cluster close, 2026-10-20 19:46) and was cancelled; time.txt is therefore empty "
            "(outFile.txt was written ~6 h 55 min after the job started at 12:50).\n")
    (out_dir / "RECOVERED.txt").write_text(note + "\n".join(lines) + "\n")
    print("\n".join(lines))
    print("GENIE3_RECOVERY_OK")


if __name__ == "__main__":
    main()
