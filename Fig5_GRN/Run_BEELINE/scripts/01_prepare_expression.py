#!/usr/bin/env python
"""Single-cell BEELINE expression input from the r9 cells of one context (common.CONTEXT: K562 or THP1).

Follows the BEELINE paper's preprocessing of experimental scRNA-seq data (Pratapa et al. 2020, Methods):
use the dataset's own normalized expression (here r9's `log1p` layer), and drop genes expressed in fewer than
10% of cells. With BEELINE_GENE_FILTER=chronopack the 10% filter is skipped and every r9 prefiltered gene is kept
(ChronoPack's own expression filter; see common.py). No ChronoPack or tradeSeq result is used.

- cells: all cells of the context's r9 prefiltered AnnData (K562: 6,857 K562 + TNF-alpha cells; THP1: the 787 THP1
  cells of the THP1 + K562 + LPS experiment, the cells of ChronoPack's r9 THP1 fit); time tags 1-12
- genes: the r9 prefiltered genes (K562 11,866; THP1 10,117), kept if counts > 0 in >= 10% of cells (beeline) or all
  kept (chronopack)
- pseudotime: the experimental time tag (1-12), as in the original metacell run

Writes output/prep/<context>[-cpfilt]/{ExpressionData.full.csv (genes x cells), PseudoTime.csv, forced_tfs.csv,
gene_filter.csv, prep_summary.json}.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402


def dense(x) -> np.ndarray:
    return x.toarray() if sp.issparse(x) else np.asarray(x)


def main() -> None:
    c.PREP_DIR.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(c.ADATA)
    if "cell_type" in adata.obs:
        print("cell types:", adata.obs["cell_type"].value_counts().to_dict())
    time = pd.to_numeric(adata.obs["time"].astype(str))
    print(f"{adata.n_obs} cells x {adata.n_vars} genes; time tags {time.value_counts().sort_index().to_dict()}")

    counts = adata.layers["counts"]
    frac = np.asarray((counts > 0).sum(axis=0)).ravel() / adata.n_obs
    keep = frac >= c.MIN_CELL_FRACTION
    gene_filter = pd.DataFrame({"gene": adata.var_names, "frac_cells_expressed": frac, "kept": keep})
    gene_filter.to_csv(c.PREP_DIR / "gene_filter.csv", index=False)

    # r9's log1p layer: report what it is (log1p of the cpm layer?) so the README can state it.
    sub = slice(0, min(200, adata.n_vars))
    log1p_cpm_match = bool(np.allclose(dense(adata.layers["log1p"][:, sub]), np.log1p(dense(adata.layers["cpm"][:, sub]))))

    forced = pd.DataFrame({"TF": c.FORCED_TFS})
    forced["in_r9_genes"] = forced["TF"].isin(adata.var_names)
    forced["frac_cells_expressed"] = forced["TF"].map(gene_filter.set_index("gene")["frac_cells_expressed"])
    forced["passes_10pct_filter"] = forced["frac_cells_expressed"].fillna(0) >= 0.10  # informational in both modes
    forced["passes_gene_filter"] = forced["TF"].map(gene_filter.set_index("gene")["kept"]).fillna(False).astype(bool)
    forced.to_csv(c.PREP_DIR / "forced_tfs.csv", index=False)
    print(forced.to_string(index=False))
    missing = forced.loc[~forced["in_r9_genes"], "TF"].tolist()
    if missing:  # no expression to infer from (THP1: ELK1 is not among the r9 prefiltered genes); recorded, skipped
        print(f"WARNING: forced TFs not in the {c.CONTEXT} r9 genes, skipped: {missing}")

    # Forced TFs are kept in the full matrix even when below the 10% filter; generateExpInputs.py never selects
    # them unless they pass its own tests, and step 03 adds them afterwards.
    rows = keep | adata.var_names.isin(c.FORCED_TFS)
    expr = pd.DataFrame(dense(adata.layers["log1p"][:, rows]).T, index=adata.var_names[rows], columns=adata.obs_names)
    expr.index.name = None
    expr.to_csv(c.PREP_DIR / "ExpressionData.full.csv")
    pd.DataFrame({"PseudoTime": time.to_numpy()}, index=pd.Index(adata.obs_names, name="")).to_csv(c.PREP_DIR / "PseudoTime.csv")

    summary = {
        "context": c.CONTEXT, "gene_filter": c.GENE_FILTER, "adata": str(c.ADATA.resolve()), "n_cells": int(adata.n_obs), "n_genes_r9": int(adata.n_vars),
        "min_cell_fraction": c.MIN_CELL_FRACTION, "n_genes_pass_filter": int(keep.sum()),
        "n_rows_written": int(rows.sum()), "forced_tfs_not_in_genes": missing,
        "forced_tfs_below_filter": forced.loc[forced["in_r9_genes"] & ~forced["passes_gene_filter"], "TF"].tolist(),
        "log1p_layer_equals_log1p_cpm": log1p_cpm_match, "uns_log1p": {k: str(v) for k, v in adata.uns.get("log1p", {}).items()},
    }
    (c.PREP_DIR / "prep_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("PREP_OK")


if __name__ == "__main__":
    main()
