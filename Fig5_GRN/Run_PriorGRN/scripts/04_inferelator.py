#!/usr/bin/env python
"""Inferelator 3 GRN with the CellOracle hg38 promoter base GRN as the prior (run with .conda_env).

TFA workflow + BBSR (Inferelator's default single-task regression; Gibbs et al. 2022) on the r9 log1p layer (already
normalized and log-transformed, so no further single-cell preprocessing), prior = output/inputs/<dataset>/prior.tsv,
TFs = the prior's TFs. No gold standard is given (use_no_gold_standard): Inferelator's precision/recall scoring is
skipped, so no Fig5 ground truth enters the run. 10 bootstraps, seed 777.
EdgeWeight = combined_confidences (Inferelator's edge ranking); every edge with confidence > 0 is kept.
Usage: python scripts/04_inferelator.py <dataset> <n_procs>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
from inferelator import MPControl, inferelator_workflow

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402


def main() -> None:
    ds, n_procs = sys.argv[1], int(sys.argv[2])
    d = c.NETWORK_DIR / ds / "INFERELATOR"
    d.mkdir(parents=True, exist_ok=True)
    inp = c.INPUTS_DIR / ds
    MPControl.set_multiprocess_engine("multiprocessing")
    MPControl.set_processes(n_procs)
    MPControl.connect()

    w = inferelator_workflow(regression="bbsr", workflow="tfa")
    w.set_file_paths(input_dir=str(inp), output_dir=str(d / "inferelator_output"),
                     tf_names_file="tfs_inferelator.txt", priors_file="prior.tsv")
    w.set_expression_file(h5ad="expr.h5ad", h5_layer="log1p")
    w.set_file_properties(expression_matrix_columns_are_genes=True)
    w.set_network_data_flags(use_no_gold_standard=True)
    w.set_run_parameters(num_bootstraps=10, random_seed=c.SEED)
    res = w.run()

    net = res.network.copy()
    net.to_csv(d / "network_all.tsv", sep="\t", index=False)
    e = net.rename(columns={"regulator": "Gene1", "target": "Gene2", "combined_confidences": "EdgeWeight"})
    e = e[e.EdgeWeight > 0]
    c.write_ranked_edges(e, d / "rankedEdges.csv")
    s = {"dataset": ds, "edges": int(len(e)), "tfs": int(e.Gene1.nunique()), "targets": int(e.Gene2.nunique()),
         "fig5_tf_edges": {t: int((e.Gene1 == t).sum()) for t in c.FIG5_TFS}}
    (d / "summary.json").write_text(json.dumps(s, indent=2))
    print(json.dumps(s))
    print("INFERELATOR_OK")


if __name__ == "__main__":
    main()
