#!/usr/bin/env python
"""BEELINE TFs + 1000 and TFs + 2000 datasets, plus one BLRunner config per (dataset, algorithm).

1. refNetwork source: the K562 ENCODE ChIP ground truth (Fig5) as a Gene1 (TF) -> Gene2 (target) edge list. BEELINE
   only uses it for its own evaluation; Fig5 evaluates the predictions itself.
2. BEELINE's generateExpInputs.py, unchanged, with its p-value filter off (-p 0): the paper's Bonferroni P < 0.01
   time test keeps only 277 genes on K562 (output/prep/K562/GeneOrdering.csv), too few for TFs + 1000 / 2000. Selected:
   every BEELINE TF (human-tfs.csv) in GeneOrdering.csv plus the N non-TF genes with the highest variance:
       generateExpInputs.py -e ExpressionData.full.csv -g GeneOrdering.csv -i human-tfs.csv -f chip_network.csv
                            -p 0 -t -n N -o <dataset>/
3. Fig5 TFs (common.FORCED_TFS) not selected in step 2 are appended to ExpressionData.csv (the only deviation),
   and refNetwork.csv is re-restricted to the final genes.
   With BEELINE_PVAL > 0 (common.py), step 2 runs with -p <P> instead (uncorrected; no -c): the paper's selection,
   significant TFs + the N most-varying significant genes (all of them if fewer than N).
4. Configs: beeline/BEELINE-data/config-files/scRNA-Seq/<dataset>/<ALGO>.yaml, same algorithm settings as the
   original metacell run (GENIE3, GRNBOOST2 defaults; PPCOR pearson, pVal 0.01; SINCERITIES nBins 10).

Runs for one context (common.CONTEXT). Writes output/prep/<context>/chip_network.csv, output/gene_selection_<dataset>.csv
and the context's rows of output/inputs_summary.csv (other contexts' rows are kept). The refNetwork source is the
K562 ChIP network for both contexts (Fig5 also scores THP1 methods against it).
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

BEELINE_PYTHON = os.environ.get("BEELINE_PYTHON", "/new-stg/home/alex/miniconda3/envs/BEELINE/bin/python")
ALGORITHMS = {
    "GENIE3": {"should_run": [True]},
    "GRNBOOST2": {"should_run": [True]},
    "PPCOR": {"pVal": [0.01], "corrType": ["pearson"], "should_run": [True]},
    "SINCERITIES": {"nBins": [10], "should_run": [True]},
}
CONFIG = """input_settings:
  input_dir: "inputs"
  dataset_dir: "{dataset_dir}"
  datasets:
    - name: "{name}"
      exprData: "ExpressionData.csv"
      cellData: "PseudoTime.csv"
      trueEdges: "refNetwork.csv"
  algorithms:
    - name: "{algo}"
      params:
{params}
output_settings:
  output_dir: "outputs"
  output_prefix: "{name}"
"""


def chip_network() -> Path:
    path = c.PREP_DIR / "chip_network.csv"
    gt = pd.read_csv(c.CHIP_GT, sep="\t", index_col=0)  # targets x TFs, nonzero = edge
    long = gt.stack()
    edges = long[long != 0].reset_index()
    edges.columns = ["Gene2", "Gene1", "value"]
    edges[["Gene1", "Gene2"]].to_csv(path, index=False)
    print(f"ChIP network: {len(edges):,} edges, {edges.Gene1.nunique()} TFs")
    return path


def main() -> None:
    net = chip_network()
    full = pd.read_csv(c.PREP_DIR / "ExpressionData.full.csv", index_col=0)
    order = pd.read_csv(c.PREP_DIR / "GeneOrdering.csv", index_col=0)
    tfs = set(pd.read_csv(c.HUMAN_TFS)["TF"])
    bf = 0.01 / len(order)  # the paper's cutoff, reported only
    sig = order.index[order["VGAMpValue"] < bf]
    rows = []
    for n in c.N_GENES:
        name = c.dataset_name(n)
        out_dir = c.BEELINE_ROOT / "BEELINE-data" / "inputs" / c.DATASET_DIR / name
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = [BEELINE_PYTHON, str(c.BEELINE_CODE / "generateExpInputs.py"),
               "-e", str(c.PREP_DIR / "ExpressionData.full.csv"), "-g", str(c.PREP_DIR / "GeneOrdering.csv"),
               "-i", str(c.HUMAN_TFS), "-f", str(net), "-p", str(c.PVAL), "-t", "-n", str(n),
               "-o", str(out_dir) + "/"]
        print(" ".join(cmd), flush=True)
        log = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
        (out_dir / "generateExpInputs.log").write_text(log)

        expr = pd.read_csv(out_dir / "ExpressionData.csv", index_col=0)
        selected = list(expr.index)
        added = [tf for tf in c.FORCED_TFS if tf not in expr.index and tf in full.index]  # absent genes: see forced_tfs.csv
        expr = pd.concat([expr, full.loc[added]])
        expr.to_csv(out_dir / "ExpressionData.csv")
        netdf = pd.read_csv(net)
        netdf = netdf[netdf.Gene1.isin(expr.index) & netdf.Gene2.isin(expr.index) & (netdf.Gene1 != netdf.Gene2)]
        netdf.drop_duplicates().to_csv(out_dir / "refNetwork.csv", index=False)

        sel = pd.DataFrame({"gene": expr.index})
        sel["is_beeline_tf"] = sel.gene.isin(tfs)
        sel["bonferroni_significant"] = sel.gene.isin(sig)
        sel["reason"] = ["forced_fig5_tf" if g in added else ("beeline_tf" if g in tfs else "top_variance") for g in sel.gene]
        if c.PVAL:  # everything BEELINE selected must pass the cut-off
            assert (sel.loc[sel.reason != "forced_fig5_tf", "gene"].map(order["VGAMpValue"]) < c.PVAL).all()
        sel = sel.merge(order[["VGAMpValue", "Variance"]], left_on="gene", right_index=True, how="left")
        sel.to_csv(c.OUTPUT_DIR / f"gene_selection_{name}.csv", index=False)
        rows.append({"dataset": name, "pval_cutoff_uncorrected": c.PVAL, "n_cells": expr.shape[1], "n_genes": expr.shape[0],
                     "n_beeline_tfs": int((sel.reason == "beeline_tf").sum()),
                     "n_top_variance": int((sel.reason == "top_variance").sum()),
                     "n_forced_fig5_tfs": len(added), "forced_added": ",".join(added),
                     "n_ref_edges": len(netdf), "n_tested_genes": len(order), "n_bonferroni_significant_genes": len(sig),
                     "bonferroni_cutoff_reported_only": bf})
        print(f"{name}: {expr.shape[0]} genes x {expr.shape[1]} cells ({len(selected)} selected by BEELINE, "
              f"{len(added)} forced: {added})")

        cfg_dir = c.BEELINE_ROOT / "BEELINE-data" / "config-files" / c.DATASET_DIR / name
        cfg_dir.mkdir(parents=True, exist_ok=True)
        for algo, params in ALGORITHMS.items():
            p = "\n".join(f"        {k}: {v}" for k, v in params.items())
            (cfg_dir / f"{algo}.yaml").write_text(CONFIG.format(dataset_dir=c.DATASET_DIR, name=name, algo=algo, params=p))
    summary_path = c.OUTPUT_DIR / "inputs_summary.csv"
    summary = pd.DataFrame(rows)
    if summary_path.exists():  # keep the other contexts' rows
        old = pd.read_csv(summary_path)
        summary = pd.concat([old[~old["dataset"].isin(summary["dataset"])], summary], ignore_index=True)
    summary.sort_values("dataset").to_csv(summary_path, index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    print("INPUTS_OK")


if __name__ == "__main__":
    main()
