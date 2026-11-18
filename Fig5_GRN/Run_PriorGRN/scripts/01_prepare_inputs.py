#!/usr/bin/env python
"""Per-dataset inputs for the motif-prior methods, from the Run_BEELINE P < 0.05 datasets (run with .conda_env).

Genes and cells: exactly those of Run_BEELINE's <dataset>/ExpressionData.csv (genes x cells, r9 log1p layer); the
matrix is checked against the r9 AnnData's log1p layer before anything is written.

Writes output/inputs/<dataset>/:
  expr.h5ad              cells x genes; X = r9 cpm layer (CellOracle imports it as "raw" counts and applies log1p,
                         giving the log1p layer), layers log1p and counts; obs time
  expr_log1p.tsv         cells x genes log1p matrix (pySCENIC GRNBoost2 and cisTarget module filtering)
  tfs_pyscenic.txt       pySCENIC's hs_hgnc_tfs.txt TFs among the dataset's genes (GRNBoost2 regulators)
  prior.tsv              genes x TFs 0/1 prior from the CellOracle hg38 promoter base GRN (a TF -> gene entry is 1 when
                         any of the gene's promoter windows has the TF's motif); TFs among the dataset's genes with at
                         least one target (Inferelator prior; CellOracle reads the base GRN itself)
  tfs_inferelator.txt    the prior's TFs
  summary.json
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
    base = pd.read_parquet(c.BASE_GRN)
    base_tfs = list(base.columns[2:])
    base_gene = base.drop(columns="peak_id").groupby("gene_short_name").max()  # gene x TF, 1 if any promoter window
    scenic_tfs = set(pd.read_csv(c.PYSCENIC_TFS, header=None)[0])
    print(f"base GRN: {base.shape[0]} promoter windows, {base_gene.shape[0]} genes, {len(base_tfs)} TFs; "
          f"pySCENIC TF list: {len(scenic_tfs)}")
    summaries = {}
    for ctx, datasets in c.DATASETS.items():
        adata = ad.read_h5ad(c.ADATAS[ctx])
        cpm_sums = np.asarray(adata.layers["cpm"].sum(axis=1)).ravel()
        print(f"{ctx}: {adata.n_obs} cells x {adata.n_vars} genes; cpm row sums {cpm_sums.min():.4g}-{cpm_sums.max():.4g}")
        for ds in datasets:
            expr = pd.read_csv(c.BEELINE_INPUTS / ds / "ExpressionData.csv", index_col=0)
            genes = list(expr.index)
            assert list(expr.columns) == list(adata.obs_names), f"{ds}: cells differ from the r9 AnnData"
            sub = adata[:, genes]
            log1p = dense(sub.layers["log1p"]).astype(np.float64)
            err = float(np.abs(log1p - expr.to_numpy().T).max())
            assert err < 1e-5, f"{ds}: ExpressionData.csv differs from the r9 log1p layer (max |diff| {err})"

            out = c.INPUTS_DIR / ds
            out.mkdir(parents=True, exist_ok=True)
            h = ad.AnnData(X=sp.csr_matrix(dense(sub.layers["cpm"]).astype(np.float32)),
                           obs=pd.DataFrame({"time": pd.to_numeric(sub.obs["time"].astype(str)).to_numpy()},
                                            index=sub.obs_names),
                           var=pd.DataFrame(index=pd.Index(genes)))
            h.layers["log1p"] = sp.csr_matrix(log1p.astype(np.float32))
            h.layers["counts"] = sp.csr_matrix(dense(sub.layers["counts"]))
            h.write_h5ad(out / "expr.h5ad")
            pd.DataFrame(log1p, index=sub.obs_names, columns=genes).to_csv(out / "expr_log1p.tsv", sep="\t",
                                                                          float_format="%.7g")

            gs = set(genes)
            tfs_scenic = sorted(gs & scenic_tfs)
            (out / "tfs_pyscenic.txt").write_text("\n".join(tfs_scenic) + "\n")
            prior = base_gene.reindex(index=genes, columns=[t for t in base_tfs if t in gs]).fillna(0).astype(int)
            prior = prior.loc[:, prior.sum(axis=0) > 0]
            prior.to_csv(out / "prior.tsv", sep="\t")
            (out / "tfs_inferelator.txt").write_text("\n".join(prior.columns) + "\n")

            s = {"context": ctx, "n_cells": int(h.n_obs), "n_genes": len(genes),
                 "max_abs_diff_vs_r9_log1p": err,
                 "n_tfs_pyscenic": len(tfs_scenic),
                 "n_tfs_base_grn_prior": int(prior.shape[1]),
                 "n_genes_in_base_grn": int(len(gs & set(base_gene.index))),
                 "prior_edges": int(prior.to_numpy().sum()),
                 "prior_density": float(prior.to_numpy().mean()),
                 "fig5_tfs_in_genes": [t for t in c.FIG5_TFS if t in gs],
                 "fig5_tfs_missing_pyscenic": [t for t in c.FIG5_TFS if t in gs and t not in tfs_scenic],
                 "fig5_tfs_missing_prior": [t for t in c.FIG5_TFS if t in gs and t not in prior.columns],
                 "fig5_tf_prior_targets": {t: int(prior[t].sum()) for t in c.FIG5_TFS if t in prior.columns}}
            (out / "summary.json").write_text(json.dumps(s, indent=2))
            summaries[ds] = s
            print(ds, json.dumps(s))
    (c.INPUTS_DIR / "summary.json").write_text(json.dumps(summaries, indent=2))
    print("PRIORGRN_PREP_OK")


if __name__ == "__main__":
    main()
